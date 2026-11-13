using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;

namespace LiTiming
{
    /// <summary>
    /// Stage C console: fidelity against the Python reference, timing on the game clock, stability and ledger fuzz.
    ///   LiTiming fidelity  DIR OUT.jsonl
    ///   LiTiming timing    DIR OUT.csv [--modes full,wet] [--only a,b] [--label L]
    ///   LiTiming stability DIR OUT.jsonl
    ///   LiTiming fuzz      FILE.lis OUT.jsonl [--seqs 1000]
    /// The game clock (owner card 2026-11-12): one water-second per 60 ticks, a recompute every 12 ticks = 0.2 water-s,
    /// stepped in CFL substeps; ms per tick is the recompute's cost over 12 ticks.
    /// </summary>
    public static class Program
    {
        const double Recompute = 0.2;       // water-seconds per 12-tick recompute
        const double DtCap = 5.0;
        static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

        public static int Main(string[] args)
        {
            if (args.Length < 2) { Console.Error.WriteLine("usage: fidelity|timing|stability|fuzz ..."); return 2; }
            Console.Error.WriteLine("runtime: " + RuntimeLabel());
            try { Process.GetCurrentProcess().PriorityClass = ProcessPriorityClass.High; } catch (Exception e) { Console.Error.WriteLine("priority: " + e.Message); }
            switch (args[0])
            {
                case "fidelity": return Fidelity(args[1], args[2]);
                case "timing": return Timing(args[1], args[2], Opt(args, "--modes", "full,wet").Split(','), Opt(args, "--only", ""), Opt(args, "--label", RuntimeLabel()));
                case "stability": return Stability(args[1], args[2], Opt(args, "--only", ""));
                case "profile":
                    {
                        var sc = SceneFile.Read(args[1]);
                        var m = LocalInertialDL.FromScene(sc); m.WetOnly = Opt(args, "--modes", "wet") == "wet";
                        for (int k = 0; k < 25; k++) Advance(m, null);
                        LocalInertialDL.Profile = true;
                        long s0 = m.Steps;
                        for (int k = 0; k < 200; k++) Advance(m, null);
                        long ns = m.Steps - s0;
                        string[] names = { "select", "faces", "limit", "move", "clear", "forcing", "refresh" };
                        for (int k = 0; k < 7; k++)
                            Console.WriteLine(names[k] + " " + (LocalInertialDL.Prof[k] * 1e6 / Stopwatch.Frequency / ns).ToString("F1", Inv) + " us/step");
                        Console.WriteLine("steps " + ns + " wet " + m.WetCells(0) + " active_blocks " + m.ActiveBlocks);
                        return 0;
                    }
                case "fuzz": return Fuzz(args[1], args[2], int.Parse(Opt(args, "--seqs", "1000")));
            }
            return 2;
        }

        static string Opt(string[] a, string k, string d) { int i = Array.IndexOf(a, k); return i >= 0 && i + 1 < a.Length ? a[i + 1] : d; }

        public static string RuntimeLabel()
        {
            var mono = Type.GetType("Mono.Runtime");
            if (mono != null)
            {
                var m = mono.GetMethod("GetDisplayName", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Static);
                return "Mono " + (m != null ? m.Invoke(null, null) : "?");
            }
            return Framework();
        }

        // separate and never inlined: Unity's Mono class library may lack RuntimeInformation, and the JIT resolves a
        // method's types when it compiles it
        [System.Runtime.CompilerServices.MethodImpl(System.Runtime.CompilerServices.MethodImplOptions.NoInlining)]
        static string Framework() => System.Runtime.InteropServices.RuntimeInformation.FrameworkDescription;

        static string F(double x) => x.ToString("R", Inv);

        // ---------------------------------------------------------------------------------------------- fidelity
        static int Fidelity(string dir, string outPath)
        {
            int fails = 0;
            using (var w = new StreamWriter(outPath))
                foreach (var path in Directory.GetFiles(dir, "*.lis").OrderBy(p => p))
                {
                    var s = SceneFile.Read(path);
                    if (s.Checkpoints.Count == 0) continue;
                    var full = LocalInertialDL.FromScene(s);
                    var wet = LocalInertialDL.FromScene(s); wet.WetOnly = true;
                    var eta0 = Eta(full);
                    var digs = new Queue<(double t, bool[] mask, double[] depth)>(s.Digs);
                    var digs2 = new Queue<(double t, bool[] mask, double[] depth)>(s.Digs);
                    bool[] pit = s.Digs.Count > 0 ? s.Digs[0].mask : null;
                    int x0 = s.Nx / 3;
                    foreach (var cp in s.Checkpoints)
                    {
                        int steps = 0;
                        foreach (var m in new[] { full, wet })
                        {
                            var q = m == full ? digs : digs2;
                            while (m.T < cp.t - 1e-9)
                            {
                                while (q.Count > 0 && m.T >= q.Peek().t - 1e-9) { var d = q.Dequeue(); m.Dig(d.mask, d.depth); }
                                double dt = Math.Min(Math.Min(m.StableDt(), DtCap), cp.t - m.T);
                                if (q.Count > 0 && q.Peek().t - m.T > 1e-9) dt = Math.Min(dt, q.Peek().t - m.T);
                                m.Step(dt);
                                if (m == full) steps++;
                            }
                        }
                        var h = full.DepthArray();
                        var hp = cp.h;
                        double l1 = 0, sp = 0; int inter = 0, union = 0;
                        for (int c = 0; c < h.Length; c++)
                        {
                            l1 += Math.Abs(h[c] - hp[c]); sp += Math.Abs(hp[c]);
                            bool a = h[c] > 1e-3, b = hp[c] > 1e-3;
                            if (a && b) inter++;
                            if (a || b) union++;
                        }
                        double l1rel = sp > 0 ? l1 / sp : 0, iou = union > 0 ? (double)inter / union : 1;
                        bool same = full.V.SequenceEqual(wet.V);
                        var row = new StringBuilder();
                        row.Append("{\"scene\":\"").Append(s.Name).Append("\",\"t\":").Append(F(cp.t))
                           .Append(",\"cs_steps_total\":").Append(full.Steps).Append(",\"py_steps_total\":").Append(cp.steps)
                           .Append(",\"l1rel\":").Append(F(l1rel)).Append(",\"iou\":").Append(F(iou))
                           .Append(",\"books_dl\":").Append(full.BooksError).Append(",\"books_dl_wet\":").Append(wet.BooksError)
                           .Append(",\"min_v\":").Append(Math.Min(full.MinVSeen, wet.MinVSeen))
                           .Append(",\"wet_only_identical\":").Append(same ? "true" : "false")
                           .Append(",\"py_volume_error_m3\":").Append(F(cp.verr));
                        bool ok = l1rel <= 0.02 && iou >= 0.98 && full.BooksError == 0 && wet.BooksError == 0 && same && full.MinVSeen >= 0;
                        if (s.Name.StartsWith("lake_at_rest"))
                        {
                            double mv = 0; var e = Eta(full);
                            for (int c = 0; c < e.Length; c++) if (full.V[c] > 0) mv = Math.Max(mv, Math.Abs(e[c] - eta0[c]));
                            row.Append(",\"max_surface_move_m\":").Append(F(mv));
                            ok &= mv < 1e-3;
                        }
                        if (pit != null)
                        {
                            double a = 0, b = 0; int n = 0;
                            for (int c = 0; c < h.Length; c++) if (pit[c]) { a += h[c]; b += hp[c]; n++; }
                            double rel = b > 0 ? Math.Abs(a - b) / b : Math.Abs(a);
                            row.Append(",\"pit_mean_cs\":").Append(F(a / n)).Append(",\"pit_mean_py\":").Append(F(b / n));
                            if (cp.t >= s.Duration - 1e-9) ok &= rel <= 0.02;
                        }
                        if (s.Name.StartsWith("dam_break"))
                        {
                            double a = 0, b = 0;
                            for (int j = 0; j < s.Ny; j++) for (int i = x0 + 2; i < s.Nx; i++) { a += h[j * s.Nx + i]; b += hp[j * s.Nx + i]; }
                            row.Append(",\"east_vol_rel_err\":").Append(F(Math.Abs(a - b) / b));
                            ok &= Math.Abs(a - b) / b <= 0.02;
                        }
                        row.Append(",\"pass\":").Append(ok ? "true" : "false").Append('}');
                        w.WriteLine(row); w.Flush();
                        Console.WriteLine(row);
                        if (!ok) fails++;
                    }
                }
            Console.WriteLine(fails == 0 ? "FIDELITY PASS" : "FIDELITY FAIL rows=" + fails);
            return fails == 0 ? 0 : 1;
        }

        static double[] Eta(LocalInertialDL m)
        {
            var e = new double[m.V.Length];
            for (int c = 0; c < e.Length; c++) e[c] = m.Z[c] + m.V[c] * m.MPerDl;
            return e;
        }

        // ------------------------------------------------------------------------------------------------ timing
        /// <summary>Advance one recompute (0.2 water-s) in CFL substeps, applying digs due. Returns substeps.</summary>
        static int Advance(LocalInertialDL m, Queue<(double t, bool[] mask, double[] depth)> digs)
        {
            while (digs != null && digs.Count > 0 && m.T >= digs.Peek().t - 1e-9) { var d = digs.Dequeue(); m.Dig(d.mask, d.depth); }
            double end = m.T + Recompute; int k = 0;
            while (m.T < end - 1e-9)
            {
                m.Step(Math.Min(Math.Min(m.StableDt(), DtCap), end - m.T));
                k++;
            }
            return k;
        }

        static int Timing(string dir, string outPath, string[] modes, string only, string label)
        {
            var onlySet = new HashSet<string>(only.Split(new[] { ',' }, StringSplitOptions.RemoveEmptyEntries));
            bool header = !File.Exists(outPath);
            using (var w = new StreamWriter(outPath, append: true))
            {
                if (header) w.WriteLine("runtime,scene,n,mode,recomputes,substeps,med_ms_recompute,p95_ms_recompute,mean_ms_recompute,max_ms_recompute,med_ms_tick,p95_ms_tick,mean_ms_tick,mean_wet_cells,mean_active_frac,books_dl");
                foreach (var path in Directory.GetFiles(dir, "*.lis").OrderBy(p => p))
                {
                    var s = SceneFile.Read(path);
                    string scene = s.Name.Substring(0, s.Name.LastIndexOf('_'));
                    if (onlySet.Count > 0 && !onlySet.Contains(scene)) continue;
                    foreach (var mode in modes)
                    {
                        // warm-up on a throwaway copy: JIT and caches, 25 recomputes
                        var warm = LocalInertialDL.FromScene(s); warm.WetOnly = mode == "wet";
                        var wq = new Queue<(double, bool[], double[])>(s.Digs);
                        for (int k = 0; k < 25; k++) Advance(warm, wq);

                        var m = LocalInertialDL.FromScene(s); m.WetOnly = mode == "wet";
                        var q = new Queue<(double, bool[], double[])>(s.Digs);
                        int nrec = (int)Math.Round(s.Duration / Recompute);
                        var ms = new double[nrec];
                        long sub = 0; double wetSum = 0, actSum = 0;
                        int nb = ((s.Nx + LocalInertialDL.B - 1) / LocalInertialDL.B) * ((s.Ny + LocalInertialDL.B - 1) / LocalInertialDL.B);
                        for (int r = 0; r < nrec; r++)
                        {
                            long t0 = Stopwatch.GetTimestamp();
                            sub += Advance(m, q);
                            long t1 = Stopwatch.GetTimestamp();
                            ms[r] = (t1 - t0) * 1000.0 / Stopwatch.Frequency;
                            if ((r & 7) == 0) { wetSum += m.WetCells(0.0); actSum += (double)m.ActiveBlocks / nb; }
                        }
                        int samples = (nrec + 7) / 8;
                        var sorted = ms.OrderBy(x => x).ToArray();
                        double med = sorted[nrec / 2], p95 = sorted[(int)Math.Min(nrec - 1, Math.Ceiling(0.95 * nrec) - 1)];
                        double mean = ms.Average(), max = sorted[nrec - 1];
                        string row = string.Join(",", label.Replace(',', ' '), scene, s.Nx.ToString(Inv), mode, nrec.ToString(Inv), sub.ToString(Inv),
                            med.ToString("F4", Inv), p95.ToString("F4", Inv), mean.ToString("F4", Inv), max.ToString("F4", Inv),
                            (med / 12).ToString("F4", Inv), (p95 / 12).ToString("F4", Inv), (mean / 12).ToString("F4", Inv),
                            (wetSum / samples).ToString("F0", Inv), (actSum / samples).ToString("F3", Inv), m.BooksError.ToString(Inv));
                        w.WriteLine(row); w.Flush();
                        Console.WriteLine(row);
                    }
                }
            }
            return 0;
        }

        // --------------------------------------------------------------------------------------------- stability
        static int Stability(string dir, string outPath, string only)
        {
            // (scene file, water-seconds, mode): lake at rest 10,000 river-s; river 100,000 recomputes; swell and dam break long
            var runs = new List<(string file, double seconds)>
            {
                ("lake_at_rest_250.lis", 10000), ("river_250.lis", 100000 * Recompute), ("swell_250.lis", 2000),
                ("dam_break_250.lis", 2000), ("pit_dig_250.lis", 2000), ("stress_half_250.lis", 2000),
            };
            int fails = 0;
            using (var w = new StreamWriter(outPath, append: true))
                foreach (var (file, seconds) in runs)
                {
                    if (only.Length > 0 && !file.StartsWith(only)) continue;
                    var s = SceneFile.Read(Path.Combine(dir, file));
                    var m = LocalInertialDL.FromScene(s); m.WetOnly = true;
                    var q = new Queue<(double, bool[], double[])>(s.Digs);
                    var eta0 = Eta(m);
                    long nrec = (long)Math.Round(seconds / Recompute), worstBooks = 0, sub = 0;
                    double maxMove = 0, vmaxM = 0; bool finite = true;
                    var sw = Stopwatch.StartNew();
                    for (long r = 0; r < nrec; r++)
                    {
                        sub += Advance(m, q);
                        long be = m.BooksError;
                        if (be != 0 && Math.Abs(be) > Math.Abs(worstBooks)) worstBooks = be;
                        if ((r % 50) == 0 || r == nrec - 1)
                        {
                            double mx = 0;
                            for (int c = 0; c < m.V.Length; c++) mx = Math.Max(mx, m.V[c] * m.MPerDl);
                            vmaxM = Math.Max(vmaxM, mx);
                            double dt = m.StableDt(); if (double.IsNaN(dt) || double.IsInfinity(dt)) finite = false;
                            if (file.StartsWith("lake_at_rest"))
                            {
                                var e = Eta(m);
                                for (int c = 0; c < e.Length; c++) if (m.V[c] > 0) maxMove = Math.Max(maxMove, Math.Abs(e[c] - eta0[c]));
                            }
                        }
                    }
                    bool ok = worstBooks == 0 && m.MinVSeen >= 0 && finite && vmaxM < 50;
                    if (file.StartsWith("lake_at_rest")) ok &= maxMove < 1e-3;
                    string row = "{\"run\":\"" + file + "\",\"water_s\":" + F(seconds) + ",\"recomputes\":" + nrec + ",\"substeps\":" + sub
                        + ",\"worst_books_dl\":" + worstBooks + ",\"final_books_dl\":" + m.BooksError + ",\"min_v\":" + m.MinVSeen
                        + ",\"max_depth_m\":" + F(vmaxM) + ",\"finite\":" + (finite ? "true" : "false")
                        + (file.StartsWith("lake_at_rest") ? ",\"max_surface_move_m\":" + F(maxMove) : "")
                        + ",\"wall_s\":" + F(sw.Elapsed.TotalSeconds) + ",\"runtime\":\"" + RuntimeLabel() + "\",\"pass\":" + (ok ? "true" : "false") + "}";
                    w.WriteLine(row); w.Flush(); Console.WriteLine(row);
                    if (!ok) fails++;
                }
            Console.WriteLine(fails == 0 ? "STABILITY PASS" : "STABILITY FAIL runs=" + fails);
            return fails == 0 ? 0 : 1;
        }

        // -------------------------------------------------------------------------------------------------- fuzz
        /// <summary>Random dig, earth-dam, wall and demolish sequences on a running scene; books and floors checked after
        /// every recompute. Half the sequences run wet-only, a quarter add two-way groundwater and evaporation.</summary>
        static int Fuzz(string file, string outPath, int seqs)
        {
            var s = SceneFile.Read(file);
            int fails = 0; long checks = 0, edits = 0, refused = 0;
            var sw = Stopwatch.StartNew();
            using (var w = new StreamWriter(outPath))
            {
                for (int k = 0; k < seqs; k++)
                {
                    var rng = new Random(1000 + k);
                    var m = LocalInertialDL.FromScene(s);
                    m.WetOnly = (k & 1) == 1;
                    if ((k & 3) == 2) { m.Evaporation = 1e-3; m.GroundwaterBothWays = true; m.SetGroundwater(0.5, 5e-3); }
                    var walls = new List<int>();
                    string fail = null;
                    int nrec = 150;    // 30 water-s
                    for (int r = 0; r < nrec && fail == null; r++)
                    {
                        if (rng.Next(3) == 0)
                        {
                            int op = rng.Next(4), cj = rng.Next(s.Ny), ci = rng.Next(s.Nx), hw = 1 + rng.Next(6);
                            edits++;
                            for (int j = Math.Max(0, cj - hw); j < Math.Min(s.Ny, cj + hw); j++)
                                for (int i = Math.Max(0, ci - hw); i < Math.Min(s.Nx, ci + hw); i++)
                                {
                                    int c = j * s.Nx + i;
                                    if (op == 0) m.Z[c] -= 0.1 + 3 * rng.NextDouble();                  // dig
                                    else if (op == 1) m.RaiseBed(c, 0.2 + 3 * rng.NextDouble());        // earth dam, water lifted
                                    else if (op == 2 && (j == cj || i == ci)) { if (m.BuildWall(c)) walls.Add(c); else refused++; }
                                }
                            if (op == 3 && walls.Count > 0)
                                for (int x = rng.Next(1, Math.Min(walls.Count, 20) + 1); x > 0 && walls.Count > 0; x--)
                                {
                                    int ix = rng.Next(walls.Count); m.Demolish(walls[ix]); walls.RemoveAt(ix);
                                }
                        }
                        Advance(m, null);
                        checks++;
                        if (m.BooksError != 0) fail = "books " + m.BooksError;
                        else if (m.MinVSeen < 0) fail = "below empty " + m.MinVSeen;
                        else if (double.IsNaN(m.StableDt())) fail = "nan";
                    }
                    if (fail != null)
                    {
                        fails++;
                        w.WriteLine("{\"seq\":" + k + ",\"fail\":\"" + fail + "\"}"); w.Flush();
                    }
                }
                string sum = "{\"file\":\"" + Path.GetFileName(file) + "\",\"seqs\":" + seqs + ",\"recompute_checks\":" + checks + ",\"edits\":" + edits
                    + ",\"walls_refused\":" + refused + ",\"fails\":" + fails + ",\"wall_s\":" + F(sw.Elapsed.TotalSeconds) + ",\"runtime\":\"" + RuntimeLabel() + "\"}";
                w.WriteLine(sum); Console.WriteLine(sum);
            }
            Console.WriteLine(fails == 0 ? "FUZZ PASS" : "FUZZ FAIL seqs=" + fails);
            return fails == 0 ? 0 : 1;
        }
    }
}
