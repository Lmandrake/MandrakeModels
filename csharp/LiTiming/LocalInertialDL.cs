using System;
using System.Collections.Generic;

namespace LiTiming
{
    /// <summary>
    /// Local-inertial shallow water (Bates 2010, de Almeida 2012) on whole decilitres: a port of MandrakeModels
    /// <c>solvers/local_inertial.py</c> plus <c>base.py</c> transfer and forcing.
    ///
    /// Water is an int per square (dL). Each interior face keeps a double discharge (m^2/s, the momentum memory) and a
    /// rounding residue in dL: the face's float transfer plus its residue is rounded to whole dL, the remainder carried to
    /// the next step, so sub-decilitre flows still happen over time. The donor limiter is integer (each donor's outflows
    /// floor-scaled to what it holds), then each face moves its whole dL once with opposite signs: the books are exact and
    /// no square goes below empty. External water enters and leaves only through booked forcing terms.
    ///
    /// WetOnly steps only 16x16 blocks holding water, dilated by one block. Water crosses at most one square per step, so a
    /// dry block two blocks from water cannot change; the result is identical to full-grid stepping (checked by the harness).
    /// </summary>
    public sealed class LocalInertialDL
    {
        public const double G = 9.81;
        public const double DlPerM3 = 10000.0;
        const double WetFace = 1e-4;           // m; Python hx > 1e-4
        public const int B = 16;

        public readonly int Ny, Nx;
        public readonly double Dx, Area, MPerDl;  // MPerDl: metres of depth per dL in one square
        public double Manning, Theta, Alpha;
        /// <summary>Face Froude cap (MandrakeModels <c>local_inertial_frcap</c>): every face discharge, after friction, is
        /// clipped to FrMax * h_f * sqrt(g h_f), h_f the Bates face depth; the clipped value is what moves and what the face
        /// remembers. 0.544 = (2/3)^1.5 is the critical-flow (broad-crested weir) discharge. 0 = no cap.</summary>
        public double FrMax;
        public double T;
        public long Steps;
        public bool WetOnly;
        public bool GroundwaterBothWays;

        public readonly double[] Z;
        public readonly int[] V;
        public readonly bool[] Open;
        readonly bool[] openX, openY;
        // qx, qy carry the closed outer faces (always 0) so neighbour reads need no edge tests: qx is ny x (nx+1), the face
        // left of square (j,i) at j*(nx+1)+i; qy is (ny+1) x nx, the face above (j,i) at j*nx+i. The per-step arrays
        // (nqx, rx, fx, openX; nqy, ry, fy, openY) hold interior faces only: x-face f=j*(nx-1)+i sits at qx[f+2j+1],
        // y-face f=j*nx+i at qy[f+nx].
        readonly double[] qx, qy, nqx, nqy, rx, ry;
        readonly int[] fx, fy, outf, lim;

        // forcing
        public readonly List<(int[] cells, double q, double acc)> Sources = new List<(int[], double, double)>();
        public readonly List<(int[] cells, double rate)> Sinks = new List<(int[], double)>();
        public readonly List<(int[] cells, int kind, double a, double b)> Levels = new List<(int[], int, double, double)>();
        public double Evaporation, Rain;
        double evapAcc, rainAcc;
        public bool HasGroundwater;
        public double GwTable, GwK;
        double[] gwAcc;

        // books, dL
        public long BookSource, BookRain, BookGround, BookEvap, BookSink, BookLevel, BookEdit;
        public long V0;
        public int MinVSeen;
        public long LimitedFaces;
        public long CappedFaces;      // faces the Froude cap clipped, summed over every step

        // blocks
        readonly int nbx, nby;
        readonly bool[] blockWet, blockActive, blockActivePrev, blockPinned;
        readonly int[] activeList;
        int nActive;
        int vmax;

        public LocalInertialDL(int ny, int nx, double dx, double[] z, double[] h, bool[] wall,
                               double manning = 0.03, double theta = 0.8, double alpha = 0.6, double frMax = 0)
        {
            Ny = ny; Nx = nx; Dx = dx; Area = dx * dx; MPerDl = 1.0 / (Area * DlPerM3);
            Manning = manning; Theta = theta; Alpha = alpha; FrMax = frMax;
            int n = ny * nx;
            Z = (double[])z.Clone();
            V = new int[n];
            Open = new bool[n];
            for (int c = 0; c < n; c++)
            {
                Open[c] = !wall[c];
                V[c] = Open[c] ? ToDl(Math.Max(h[c], 0.0)) : 0;
            }
            openX = new bool[ny * (nx - 1)];
            openY = new bool[(ny - 1) * nx];
            qx = new double[ny * (nx + 1)]; nqx = new double[openX.Length]; rx = new double[openX.Length]; fx = new int[openX.Length];
            qy = new double[(ny + 1) * nx]; nqy = new double[openY.Length]; ry = new double[openY.Length]; fy = new int[openY.Length];
            outf = new int[n]; lim = new int[n];
            nbx = (nx + B - 1) / B; nby = (ny + B - 1) / B;
            blockWet = new bool[nbx * nby]; blockActive = new bool[nbx * nby]; blockActivePrev = new bool[nbx * nby];
            blockPinned = new bool[nbx * nby];
            activeList = new int[nbx * nby];
            RebuildOpenFaces();
            for (int b = 0; b < blockActivePrev.Length; b++) blockActivePrev[b] = true;
            V0 = TotalVolume();
            MinVSeen = int.MaxValue;
            RefreshWetAndMax(all: true);
        }

        public int ToDl(double depthM) => (int)Math.Floor(depthM * Area * DlPerM3 + 0.5);
        public double Depth(int c) => V[c] * MPerDl;

        public static LocalInertialDL FromScene(SceneFile s)
        {
            var m = new LocalInertialDL(s.Ny, s.Nx, s.Dx, s.Z, s.H, s.Wall, s.Manning, s.Theta, s.Alpha, s.FrMax);
            // planted faults, for proving the fidelity checker goes red (never set in a measured run)
            var mut = Environment.GetEnvironmentVariable("LITIMING_MUTATE") ?? "";
            if (mut == "manning") m.Manning *= 1.3;
            if (mut == "theta") m.Theta = 1.0;
            if (mut == "frcap1" && m.FrMax > 0) m.FrMax = 1.0;   // a weaker cap than the reference's 0.544
            if (mut == "nocap") m.FrMax = 0;                      // the uncapped scheme against a capped reference
            if (mut == "nocbrt") NoCbrt = true;        // timing probe only: wrong friction
            if (mut == "walldelete") MutateWallDeletes = true;
            foreach (var (mask, q) in s.Sources) m.Sources.Add((Cells(mask, m.Open), q, 0.0));
            foreach (var (mask, rate) in s.Sinks) m.Sinks.Add((Cells(mask, m.Open), rate));
            foreach (var (mask, kind, a, b) in s.Levels) m.Levels.Add((Cells(mask, m.Open), kind, a, b));
            m.Evaporation = s.Evaporation; m.Rain = s.Rain;
            if (s.HasGroundwater) m.SetGroundwater(s.GwTable, s.GwK);
            m.PinForcingBlocks();
            return m;
        }

        public static int[] Cells(bool[] mask, bool[] open)
        {
            var l = new List<int>();
            for (int c = 0; c < mask.Length; c++) if (mask[c] && open[c]) l.Add(c);
            return l.ToArray();
        }

        public void SetGroundwater(double table, double k)
        {
            HasGroundwater = true; GwTable = table; GwK = k; gwAcc = new double[Ny * Nx];
            PinForcingBlocks();
        }

        /// <summary>Blocks that forcing can wet from nothing (sources, levels, groundwater below its table, rain) stay active.</summary>
        public void PinForcingBlocks()
        {
            Array.Clear(blockPinned, 0, blockPinned.Length);
            foreach (var s in Sources) foreach (int c in s.cells) blockPinned[BlockOf(c)] = true;
            foreach (var l in Levels) foreach (int c in l.cells) blockPinned[BlockOf(c)] = true;
            if (Rain != 0) for (int b = 0; b < blockPinned.Length; b++) blockPinned[b] = true;
            if (HasGroundwater)
                for (int c = 0; c < Z.Length; c++) if (Open[c] && Z[c] < GwTable) blockPinned[BlockOf(c)] = true;
        }

        int BlockOf(int c) => (c / Nx / B) * nbx + (c % Nx) / B;

        void RebuildOpenFaces()
        {
            for (int j = 0; j < Ny; j++)
                for (int i = 0; i < Nx - 1; i++) openX[j * (Nx - 1) + i] = Open[j * Nx + i] && Open[j * Nx + i + 1];
            for (int j = 0; j < Ny - 1; j++)
                for (int i = 0; i < Nx; i++) openY[j * Nx + i] = Open[j * Nx + i] && Open[(j + 1) * Nx + i];
        }

        public long TotalVolume()
        {
            long s = 0;
            for (int c = 0; c < V.Length; c++) s += V[c];
            return s;
        }

        public long BooksNet => BookSource + BookRain + BookGround + BookEvap + BookSink + BookLevel + BookEdit;
        /// <summary>Volume minus start minus every booked term, dL. Must be exactly 0.</summary>
        public long BooksError => TotalVolume() - V0 - BooksNet;

        public double StableDt() => Alpha * Dx / Math.Sqrt(G * Math.Max(vmax * MPerDl, 1e-3));

        // ------------------------------------------------------------------------------------------------- the step
        /// <summary>Per-phase Stopwatch ticks (select, faces, limit, move, clear, forcing, refresh) when Profile is set.</summary>
        public static bool Profile, NoCbrt, MutateWallDeletes;
        public static readonly long[] Prof = new long[7];
        static long Mark(int k, long t0) { long t = System.Diagnostics.Stopwatch.GetTimestamp(); Prof[k] += t - t0; return t; }

        public void Step(double dt)
        {
            long pt = Profile ? System.Diagnostics.Stopwatch.GetTimestamp() : 0;
            SelectBlocks();
            if (Profile) pt = Mark(0, pt);
            double inv = MPerDl, n2 = Manning * Manning, th = Theta, toDl = Dx * dt * DlPerM3;
            double gdt = G * dt, gdtdx = G * dt / Dx, half1mth = 0.5 * (1 - th);
            int nx = Nx, ny = Ny, nfx = Nx - 1;
            // Math.Max and Math.Floor are avoided in the face loop: under Unity's Mono they are not inlined.
            // every array in a local: Mono reloads a field on each access and cannot hoist it out of the loop
            var z = Z; var v = V; var qx = this.qx; var qy = this.qy; var nqx = this.nqx; var nqy = this.nqy;
            var rx = this.rx; var ry = this.ry; var fx = this.fx; var fy = this.fy; var outf = this.outf; var lim = this.lim;
            var openX = this.openX; var openY = this.openY; var activeList = this.activeList; int nbx = this.nbx, nActive = this.nActive;

            // pass 1: new discharge on every face of an active block (reads only previous q), integer flux, donor outflow
            Faces(inv, n2, th, toDl, gdt, gdtdx, half1mth);
            if (Profile) pt = Mark(1, pt);
            // pass 2: donors asked for more than they hold are limited to what they hold (snapshot before any move)
            for (int a = 0; a < nActive; a++)
            {
                int b = activeList[a], bj = b / nbx, bi = b % nbx;
                int j0 = bj * B, j1 = Math.Min(j0 + B, ny), i0 = bi * B, i1 = Math.Min(i0 + B, nx);
                for (int j = j0; j < j1; j++)
                    for (int c = j * nx + i0, ce = j * nx + i1; c < ce; c++)
                        lim[c] = outf[c] > v[c] ? v[c] : -1;
            }

            if (Profile) pt = Mark(2, pt);
            // pass 3: move whole dL once per face, opposite signs; store discharge memory
            long limited = 0;
            for (int a = 0; a < nActive; a++)
            {
                int b = activeList[a], bj = b / nbx, bi = b % nbx;
                int j0 = bj * B, j1 = Math.Min(j0 + B, ny), i0 = bi * B, i1 = Math.Min(i0 + B, nx);
                for (int j = j0; j < j1; j++)
                {
                    int row = j * nx, frow = j * nfx;
                    for (int i = i0; i < i1; i++)
                    {
                        int c = row + i;
                        if (i < nfx)
                        {
                            int f = frow + i, fi = fx[f];
                            double nq = nqx[f];
                            if (fi != 0)
                            {
                                int d = fi > 0 ? c : c + 1;
                                if (lim[d] >= 0)
                                {
                                    int mag = fi > 0 ? fi : -fi;
                                    int cut = (int)((long)mag * lim[d] / outf[d]);
                                    fi = fi > 0 ? cut : -cut;
                                    // memory scaled with the cut, as the float limiter does; never fi/(dx dt), which a
                                    // whole dL rounded up from a short step turns into a huge discharge
                                    nq *= (double)cut / mag; rx[f] = 0; limited++;
                                }
                                v[c] -= fi; v[c + 1] += fi;
                            }
                            qx[f + 2 * j + 1] = nq;
                        }
                        if (j < ny - 1)
                        {
                            int f = c, fi = fy[f];
                            double nq = nqy[f];
                            if (fi != 0)
                            {
                                int d = fi > 0 ? c : c + nx;
                                if (lim[d] >= 0)
                                {
                                    int mag = fi > 0 ? fi : -fi;
                                    int cut = (int)((long)mag * lim[d] / outf[d]);
                                    fi = fi > 0 ? cut : -cut;
                                    // memory scaled with the cut, as the float limiter does; never fi/(dx dt), which a
                                    // whole dL rounded up from a short step turns into a huge discharge
                                    nq *= (double)cut / mag; ry[f] = 0; limited++;
                                }
                                v[c] -= fi; v[c + nx] += fi;
                            }
                            qy[f + nx] = nq;
                        }
                    }
                }
            }
            LimitedFaces = limited;
            if (Profile) pt = Mark(3, pt);
            for (int a = 0; a < nActive; a++)
            {
                int b = activeList[a], bj = b / nbx, bi = b % nbx;
                int j0 = bj * B, j1 = Math.Min(j0 + B, ny), i0 = bi * B, i1 = Math.Min(i0 + B, nx);
                for (int j = j0; j < j1; j++) Array.Clear(outf, j * nx + i0, i1 - i0);
            }

            if (Profile) pt = Mark(4, pt);
            ApplyForcing(dt);
            if (Profile) pt = Mark(5, pt);
            T += dt;
            Steps++;
            RefreshWetAndMax(all: false);
            if (Profile) Mark(6, pt);
        }

        /// <summary>Pass 1 through raw pointers on the padded q layout: no bounds checks and no edge tests inside a face.
        /// Unity's Mono allocates registers per basic block, so every branch removed from this loop is a spill saved.</summary>
        unsafe void Faces(double inv, double n2, double th, double toDl, double gdt, double gdtdx, double half1mth)
        {
            int nx = Nx, ny = Ny, nfx = Nx - 1, sx = Nx + 1;
            bool noCbrt = NoCbrt, cap = FrMax > 0;
            double frSqrtG = FrMax * Math.Sqrt(G);
            long capped = 0;
            fixed (double* z = Z, qx = this.qx, qy = this.qy, nqx = this.nqx, nqy = this.nqy, rx = this.rx, ry = this.ry)
            fixed (int* v = V, fx = this.fx, fy = this.fy, outf = this.outf)
            fixed (bool* openX = this.openX, openY = this.openY)
            {
                for (int a = 0; a < nActive; a++)
                {
                    int b = activeList[a], bj = b / nbx, bi = b % nbx;
                    int j0 = bj * B, j1 = j0 + B < ny ? j0 + B : ny, i0 = bi * B, i1 = i0 + B < nx ? i0 + B : nx;
                    for (int j = j0; j < j1; j++)
                    {
                        int row = j * nx, frow = j * nfx, prow = j * sx;
                        bool down = j < ny - 1;
                        for (int i = i0; i < i1; i++)
                        {
                            int c = row + i;
                            double zc = z[c], ec = zc + v[c] * inv;
                            if (i < nfx)
                            {
                                int f = frow + i;
                                if (!openX[f]) { fx[f] = 0; nqx[f] = 0; }
                                else
                                {
                                    double zr = z[c + 1], er = zr + v[c + 1] * inv;
                                    double hf = (ec > er ? ec : er) - (zc > zr ? zc : zr);
                                    if (hf <= WetFace) { fx[f] = 0; nqx[f] = 0; rx[f] = 0; }
                                    else
                                    {
                                        double* qp = qx + prow + i + 1;              // this face, padded
                                        double q0 = *qp;
                                        double qt = th * q0 + half1mth * (qp[-1] + qp[1]);
                                        // cross discharge on the face: mean of cell-centred qy at c and c+1
                                        double qc = 0.25 * (qy[c] + qy[c + nx] + qy[c + 1] + qy[c + 1 + nx]);
                                        double mag = Math.Sqrt(q0 * q0 + qc * qc);
                                        double h73 = noCbrt ? hf * hf * hf : hf * hf * Cbrt(hf);
                                        double nq = (qt - gdtdx * hf * (er - ec)) / (1 + gdt * n2 * mag / h73);
                                        if (cap)
                                        {
                                            double qm = frSqrtG * hf * Math.Sqrt(hf);
                                            if (nq > qm) { nq = qm; capped++; } else if (nq < -qm) { nq = -qm; capped++; }
                                        }
                                        nqx[f] = nq;
                                        double tt = nq * toDl + rx[f];
                                        int fi = (int)(tt >= 0 ? tt + 0.5 : tt - 0.5);
                                        rx[f] = tt - fi;
                                        fx[f] = fi;
                                        if (fi > 0) outf[c] += fi; else if (fi < 0) outf[c + 1] -= fi;
                                    }
                                }
                            }
                            if (down)
                            {
                                int f = c;
                                if (!openY[f]) { fy[f] = 0; nqy[f] = 0; }
                                else
                                {
                                    int cd = c + nx;
                                    double zd = z[cd], ed = zd + v[cd] * inv;
                                    double hf = (ec > ed ? ec : ed) - (zc > zd ? zc : zd);
                                    if (hf <= WetFace) { fy[f] = 0; nqy[f] = 0; ry[f] = 0; }
                                    else
                                    {
                                        double* qp = qy + cd;                         // this face, padded
                                        double q0 = *qp;
                                        double qt = th * q0 + half1mth * (qp[-nx] + qp[nx]);
                                        double* xp = qx + prow + i;                   // left face of c; the row below is +sx
                                        double qc = 0.25 * (xp[0] + xp[1] + xp[sx] + xp[sx + 1]);
                                        double mag = Math.Sqrt(q0 * q0 + qc * qc);
                                        double h73 = noCbrt ? hf * hf * hf : hf * hf * Cbrt(hf);
                                        double nq = (qt - gdtdx * hf * (ed - ec)) / (1 + gdt * n2 * mag / h73);
                                        if (cap)
                                        {
                                            double qm = frSqrtG * hf * Math.Sqrt(hf);
                                            if (nq > qm) { nq = qm; capped++; } else if (nq < -qm) { nq = -qm; capped++; }
                                        }
                                        nqy[f] = nq;
                                        double tt = nq * toDl + ry[f];
                                        int fi = (int)(tt >= 0 ? tt + 0.5 : tt - 0.5);
                                        ry[f] = tt - fi;
                                        fy[f] = fi;
                                        if (fi > 0) outf[c] += fi; else if (fi < 0) outf[cd] -= fi;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            CappedFaces += capped;
        }

        /// <summary>Pick the blocks to step: all (full grid) or wet blocks plus pinned ones, dilated by one block.
        /// A block leaving the active set has its face memory cleared (its faces are dry, where a full step writes 0).</summary>
        void SelectBlocks()
        {
            nActive = 0;
            for (int bj = 0; bj < nby; bj++)
                for (int bi = 0; bi < nbx; bi++)
                {
                    int b = bj * nbx + bi;
                    bool act = !WetOnly;
                    if (!act)
                        for (int dj = -1; dj <= 1 && !act; dj++)
                            for (int di = -1; di <= 1; di++)
                            {
                                int y = bj + dj, x = bi + di;
                                if (y < 0 || y >= nby || x < 0 || x >= nbx) continue;
                                int k = y * nbx + x;
                                if (blockWet[k] || blockPinned[k]) { act = true; break; }
                            }
                    if (!act && blockActivePrev[b]) ClearBlockFaces(b);
                    blockActive[b] = act;
                    blockActivePrev[b] = act;
                    if (act) activeList[nActive++] = b;
                }
        }

        void ClearBlockFaces(int b)
        {
            int bj = b / nbx, bi = b % nbx;
            int j0 = bj * B, j1 = Math.Min(j0 + B, Ny), i0 = bi * B, i1 = Math.Min(i0 + B, Nx);
            for (int j = j0; j < j1; j++)
                for (int i = i0; i < i1; i++)
                {
                    if (i < Nx - 1) { int f = j * (Nx - 1) + i; qx[f + 2 * j + 1] = 0; rx[f] = 0; }
                    if (j < Ny - 1) { int f = j * Nx + i; qy[f + Nx] = 0; ry[f] = 0; }
                }
        }

        public int ActiveBlocks => nActive;

        void RefreshWetAndMax(bool all)
        {
            int mx = 0, mn = MinVSeen;
            for (int b = 0; b < blockWet.Length; b++)
            {
                if (!all && !blockActive[b]) continue;     // an inactive block is dry and stayed dry
                int bj = b / nbx, bi = b % nbx;
                int j0 = bj * B, j1 = Math.Min(j0 + B, Ny), i0 = bi * B, i1 = Math.Min(i0 + B, Nx);
                bool wet = false;
                for (int j = j0; j < j1; j++)
                    for (int c = j * Nx + i0, ce = j * Nx + i1; c < ce; c++)
                    {
                        int x = V[c];
                        if (x > 0) { wet = true; if (x > mx) mx = x; }
                        if (x < mn) mn = x;
                    }
                blockWet[b] = wet;
            }
            vmax = mx; MinVSeen = mn;
        }

        // ------------------------------------------------------------------------------------------------ forcing
        void ApplyForcing(double dt)
        {
            double t = T + dt, dlPerM = Area * DlPerM3;
            for (int s = 0; s < Sources.Count; s++)
            {
                // a source spreads over its squares that are open now (a wall built on one takes no water)
                var (cells, q, acc) = Sources[s];
                int n = 0;
                foreach (int c in cells) if (Open[c]) n++;
                if (n == 0) continue;
                double want = q * dt * DlPerM3 + acc;
                long tot = (long)Math.Floor(want);
                Sources[s] = (cells, q, want - tot);
                long each = tot / n, rem = tot % n, start = Steps % n, k = 0;
                foreach (int c in cells)
                {
                    if (!Open[c]) continue;
                    long r = (k - start + n) % n;
                    V[c] += (int)(each + (r < rem ? 1 : 0));
                    k++;
                }
                BookSource += tot;
            }
            if (Rain != 0)
            {
                double want = Rain * dt * dlPerM + rainAcc;
                int each = (int)Math.Floor(want);
                rainAcc = want - each;
                if (each > 0)
                {
                    long n = 0;
                    for (int c = 0; c < V.Length; c++) if (Open[c]) { V[c] += each; n++; }
                    BookRain += n * each;
                }
            }
            if (HasGroundwater)
            {
                long tot = 0;
                for (int c = 0; c < V.Length; c++)
                {
                    if (!Open[c]) continue;
                    double eta = Z[c] + V[c] * MPerDl, d = GwTable - eta;
                    if (d < 0 && !GroundwaterBothWays) { gwAcc[c] = 0; continue; }
                    double add = GwK * d * dt;
                    if (d > 0) add = Math.Min(add, d); else add = Math.Max(add, d);
                    double want = add * dlPerM + gwAcc[c];
                    int k = (int)Math.Floor(want + 0.5);
                    if (k < -V[c]) k = -V[c];
                    gwAcc[c] = Math.Max(-1, Math.Min(1, want - k));
                    V[c] += k; tot += k;
                }
                BookGround += tot;
            }
            if (Evaporation != 0)
            {
                double want = Evaporation * dt * dlPerM + evapAcc;
                int each = (int)Math.Floor(want);
                evapAcc = want - each;
                if (each > 0)
                {
                    long tot = 0;
                    for (int a = 0; a < nActive; a++)
                    {
                        int b = activeList[a], bj = b / nbx, bi = b % nbx;
                        int j0 = bj * B, j1 = Math.Min(j0 + B, Ny), i0 = bi * B, i1 = Math.Min(i0 + B, Nx);
                        for (int j = j0; j < j1; j++)
                            for (int c = j * Nx + i0, ce = j * Nx + i1; c < ce; c++)
                            {
                                int take = V[c] < each ? V[c] : each;
                                V[c] -= take; tot += take;
                            }
                    }
                    BookEvap -= tot;
                }
            }
            foreach (var (cells, rate) in Sinks)
            {
                long tot = 0;
                bool all = double.IsInfinity(rate);
                int cap = all ? int.MaxValue : (int)Math.Floor(rate * dt * dlPerM);
                foreach (int c in cells) { if (!Open[c]) continue; int take = Math.Min(V[c], cap); V[c] -= take; tot += take; }
                BookSink -= tot;
            }
            foreach (var (cells, kind, a, b) in Levels)
            {
                double level = kind == 0 ? a : a * Math.Sin(2 * Math.PI * t / b);
                long tot = 0;
                foreach (int c in cells)
                {
                    if (!Open[c]) continue;
                    int target = ToDl(Math.Max(level - Z[c], 0.0));
                    tot += target - V[c];
                    V[c] = target;
                }
                BookLevel += tot;
            }
        }

        // -------------------------------------------------------------------------------------------- terrain edits
        /// <summary>Lower the bed; each square keeps its water (its surface drops).</summary>
        public void Dig(bool[] mask, double[] depth)
        {
            for (int c = 0; c < Z.Length; c++) if (mask[c]) Z[c] -= depth[c];
            RefreshWetAndMax(all: true);
        }

        public void RaiseBed(int c, double height) { Z[c] += height; }

        /// <summary>Wall a square. Its water is pushed to open neighbours (nearest first, breadth-first), never deleted;
        /// if no open square is reachable the wall is refused.</summary>
        public bool BuildWall(int c0)
        {
            if (!Open[c0]) return false;
            int w = V[c0];
            Open[c0] = false;
            if (MutateWallDeletes) { V[c0] = 0; w = 0; }     // planted fault for the fuzz checker only
            if (w > 0)
            {
                var seen = new HashSet<int> { c0 };
                var q = new Queue<int>(); q.Enqueue(c0);
                var ring = new List<int>();
                while (q.Count > 0 && ring.Count == 0)
                {
                    int n = q.Count;
                    for (int k = 0; k < n; k++)
                    {
                        int c = q.Dequeue(), j = c / Nx, i = c % Nx;
                        foreach (int nb in new[] { i > 0 ? c - 1 : -1, i < Nx - 1 ? c + 1 : -1, j > 0 ? c - Nx : -1, j < Ny - 1 ? c + Nx : -1 })
                        {
                            if (nb < 0 || !seen.Add(nb)) continue;
                            if (Open[nb]) ring.Add(nb); else q.Enqueue(nb);
                        }
                    }
                }
                if (ring.Count == 0) { Open[c0] = true; return false; }
                int each = w / ring.Count, rem = w % ring.Count;
                for (int k = 0; k < ring.Count; k++) V[ring[k]] += each + (k < rem ? 1 : 0);
                V[c0] = 0;
            }
            AfterWallChange(c0);
            return true;
        }

        /// <summary>Remove a wall: the square opens (empty: a wall never holds water).</summary>
        public void Demolish(int c)
        {
            if (Open[c]) return;
            Open[c] = true;
            AfterWallChange(c);
        }

        void AfterWallChange(int c)
        {
            int j = c / Nx, i = c % Nx;
            if (i < Nx - 1) { int f = j * (Nx - 1) + i; openX[f] = Open[c] && Open[c + 1]; qx[f + 2 * j + 1] = 0; rx[f] = 0; }
            if (i > 0) { int f = j * (Nx - 1) + i - 1; openX[f] = Open[c - 1] && Open[c]; qx[f + 2 * j + 1] = 0; rx[f] = 0; }
            if (j < Ny - 1) { int f = c; openY[f] = Open[c] && Open[c + Nx]; qy[f + Nx] = 0; ry[f] = 0; }
            if (j > 0) { int f = c - Nx; openY[f] = Open[c - Nx] && Open[c]; qy[f + Nx] = 0; ry[f] = 0; }
            RefreshWetAndMax(all: true);
        }

        // ------------------------------------------------------------------------------------------------- helpers
        /// <summary>Cube root for x &gt; 0: bit-trick start (~3%) and two Halley steps (cubic convergence, ~1e-14 relative).
        /// Math.Cbrt is not in .NET Framework 4.7.2 or Unity's Mono profile; Math.Pow and Newton chains cost several times more.</summary>
        static double Cbrt(double x)
        {
            double y = BitConverter.Int64BitsToDouble(BitConverter.DoubleToInt64Bits(x) / 3 + 0x2A9F7893782DA1CEL);
            double y3 = y * y * y;
            y *= (y3 + 2 * x) / (2 * y3 + x);
            y3 = y * y * y;
            return y * (y3 + 2 * x) / (2 * y3 + x);
        }

        public double[] DepthArray()
        {
            var h = new double[V.Length];
            for (int c = 0; c < V.Length; c++) h[c] = V[c] * MPerDl;
            return h;
        }

        public int WetCells(double thresholdM)
        {
            int n = 0, k = (int)Math.Ceiling(thresholdM / MPerDl);
            for (int c = 0; c < V.Length; c++) if (V[c] > k) n++;
            return n;
        }
    }
}
