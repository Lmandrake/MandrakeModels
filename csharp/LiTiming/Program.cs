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
    ///   LiTiming fidelity  DIR OUT.jsonl     (the scene files carry the reference solver's face Froude cap; LIS1 = none)
    ///   LiTiming timing    DIR OUT.csv [--modes full,wet] [--only a,b] [--label L]
    ///   LiTiming stability DIR OUT.jsonl
    ///   LiTiming fuzz      FILE.lis OUT.jsonl [--seqs 1000]
    ///   LiTiming spread    DIR OUT.csv [--modes wet,sleep] [--only a,b]   per-tick cost with each recompute spread over
    ///                      its 12 ticks, checked bitwise against the recompute in one go
    ///   LiTiming settled   DIR OUT.csv [--only lake_at_rest]  sleeping: recomputes until asleep, then cost per tick
    ///   LiTiming wake      FILE.lis OUT.jsonl  a dig into sleeping water must flow as it does in water that never slept
    ///   LiTiming evensteps FILE.lis          every driver takes equal steps (the stepping rule, LocalInertialDL.EvenDt)
    /// Modes: full (every square), wet (wet blocks, dilated), sleep (wet, and settled water sleeps).
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
                case "stability": return Stability(args[1], args[2], Opt(args, "--only", ""), Opt(args, "--modes", "wet"));
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
                case "drift":
                    {
                        // lake-at-rest drift over time: max surface movement every 100 water-s
                        var sc = SceneFile.Read(args[1]);
                        double T = double.Parse(Opt(args, "--seconds", "10000"), Inv);
                        var m = LocalInertialDL.FromScene(sc); m.WetOnly = Opt(args, "--modes", "wet") == "wet";
                        bool pydt = Array.IndexOf(args, "--pydt") >= 0;
                        var eta0 = Eta(m); double next = double.Parse(Opt(args, "--every", "100"), Inv), every = next;
                        while (m.T < T - 1e-9)
                        {
                            if (pydt) StepTo(m, null, next);
                            else Advance(m, null);
                            if (m.T >= next - 1e-9)
                            {
                                var e = Eta(m); double mv = 0; int at = -1;
                                for (int c = 0; c < e.Length; c++) if (m.V[c] > 0 && Math.Abs(e[c] - eta0[c]) > mv) { mv = Math.Abs(e[c] - eta0[c]); at = c; }
                                string where = at < 0 ? "" : " at j=" + at / m.Nx + " i=" + at % m.Nx + " z=" + m.Z[at].ToString("F3", Inv) + " v=" + m.V[at] + " eta0=" + eta0[at].ToString("F4", Inv);
                                Console.WriteLine("t=" + m.T.ToString("F1", Inv) + " max_move_m=" + mv.ToString("E3", Inv) + where + " books=" + m.BooksError + " limited=" + m.LimitedFaces);
                                next += every;
                            }
                        }
                        return 0;
                    }
                case "fuzz": return Fuzz(args[1], args[2], int.Parse(Opt(args, "--seqs", "1000")));
                case "spread": return SpreadTiming(args[1], args[2], Opt(args, "--modes", "wet,sleep").Split(','), Opt(args, "--only", ""));
                case "settled": return Settled(args[1], args[2], Opt(args, "--only", "lake_at_rest"));
                case "wake": return Wake(args[1], args[2]);
                case "evensteps": return EvenSteps(args[1]);
            }
            return 2;
        }

        static LocalInertialDL Make(SceneFile s, string mode)
        {
            var m = LocalInertialDL.FromScene(s);
            m.WetOnly = mode != "full";
            m.Sleep = mode == "sleep";
            return m;
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
                    var wet = Make(s, "wet");
                    var slp = Make(s, "sleep");
                    var eta0 = Eta(full);
                    var digs = new Queue<(double t, bool[] mask, double[] depth)>(s.Digs);
                    var digs2 = new Queue<(double t, bool[] mask, double[] depth)>(s.Digs);
                    var digs3 = new Queue<(double t, bool[] mask, double[] depth)>(s.Digs);
                    bool[] pit = s.Digs.Count > 0 ? s.Digs[0].mask : null;
                    int x0 = s.Nx / 3;
                    foreach (var cp in s.Checkpoints)
                    {
                        int steps = 0;
                        foreach (var m in new[] { full, wet, slp })
                        {
                            var q = m == full ? digs : m == wet ? digs2 : digs3;
                            int n = StepTo(m, q, cp.t);
                            if (m == full) steps += n;
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
                        var (l1s, ious) = Compare(slp.DepthArray(), hp);
                        bool same = full.V.SequenceEqual(wet.V);
                        var row = new StringBuilder();
                        row.Append("{\"scene\":\"").Append(s.Name).Append("\",\"t\":").Append(F(cp.t))
                           .Append(",\"cs_steps_total\":").Append(full.Steps).Append(",\"py_steps_total\":").Append(cp.steps)
                           .Append(",\"l1rel\":").Append(F(l1rel)).Append(",\"iou\":").Append(F(iou))
                           .Append(",\"books_dl\":").Append(full.BooksError).Append(",\"books_dl_wet\":").Append(wet.BooksError)
                           .Append(",\"min_v\":").Append(Math.Min(full.MinVSeen, wet.MinVSeen))
                           .Append(",\"wet_only_identical\":").Append(same ? "true" : "false")
                           .Append(",\"py_volume_error_m3\":").Append(F(cp.verr))
                           .Append(",\"fr_max\":").Append(F(full.FrMax)).Append(",\"capped_faces\":").Append(full.CappedFaces)
                           .Append(",\"l1rel_sleep\":").Append(F(l1s)).Append(",\"iou_sleep\":").Append(F(ious))
                           .Append(",\"books_dl_sleep\":").Append(slp.BooksError).Append(",\"sleeping_blocks\":").Append(slp.SleepingBlocks);
                        bool ok = l1rel <= 0.02 && iou >= 0.98 && full.BooksError == 0 && wet.BooksError == 0 && same && full.MinVSeen >= 0;
                        // the sleeping twin is not bitwise (settled water stops being stepped) but must meet the same tolerances
                        ok &= l1s <= 0.02 && ious >= 0.98 && slp.BooksError == 0 && slp.MinVSeen >= 0;
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

        static (double l1rel, double iou) Compare(double[] h, double[] hp)
        {
            double l1 = 0, sp = 0; int inter = 0, union = 0;
            for (int c = 0; c < h.Length; c++)
            {
                l1 += Math.Abs(h[c] - hp[c]); sp += Math.Abs(hp[c]);
                bool a = h[c] > 1e-3, b = hp[c] > 1e-3;
                if (a && b) inter++;
                if (a || b) union++;
            }
            return (sp > 0 ? l1 / sp : 0, union > 0 ? (double)inter / union : 1);
        }

        static double[] Eta(LocalInertialDL m)
        {
            var e = new double[m.V.Length];
            for (int c = 0; c < e.Length; c++) e[c] = m.Z[c] + m.V[c] * m.MPerDl;
            return e;
        }

        /// <summary>Advance to time t (a checkpoint or sample) by the stepping rule, applying digs due on the way: a dig
        /// ends one interval and starts the next. Returns steps.</summary>
        static int StepTo(LocalInertialDL m, Queue<(double t, bool[] mask, double[] depth)> q, double t)
        {
            int n = 0;
            while (m.T < t - 1e-9)
            {
                while (q != null && q.Count > 0 && m.T >= q.Peek().t - 1e-9) { var d = q.Dequeue(); m.Dig(d.mask, d.depth); }
                double end = q != null && q.Count > 0 ? Math.Min(t, q.Peek().t) : t;
                n += m.AdvanceTo(end, DtCap);
            }
            return n;
        }

        // -------------------------------------------------------------------------------------------- even steps
        /// <summary>Every driver of the engine takes equal steps (the stepping rule: LocalInertialDL.EvenDt). Runs each
        /// driver on the scene for 12 water-s, the checkpoint driver on a 0.3 s frame clock (the clock that striped the
        /// compare pages), and allows each at most one clipped tail step per dig in the window (an edit ends its interval).
        /// A ragged feed (full steps, then a tail clipped to each frame) must read above that, or the audit is blind.</summary>
        static int EvenSteps(string scenePath)
        {
            var s = SceneFile.Read(scenePath);
            const double T = 12.0, Frame = 0.3;
            int frames = (int)Math.Round(T / Frame), recs = (int)Math.Round(T / Recompute);
            var rows = new List<(string name, LocalInertialDL m, bool wantZero)>();
            var bad = Make(s, "wet");
            for (int k = 1; k <= frames; k++)
                while (bad.T < k * Frame - 1e-9) bad.Step(Math.Min(Math.Min(bad.StableDt(), DtCap), k * Frame - bad.T));
            rows.Add(("ragged_feed_control", bad, false));
            var a = Make(s, "wet"); var qa = new Queue<(double, bool[], double[])>(s.Digs);
            for (int k = 1; k <= frames; k++) StepTo(a, qa, k * Frame);
            rows.Add(("checkpoint_driver(fidelity,drift)", a, true));
            var b = Make(s, "wet"); var qb = new Queue<(double, bool[], double[])>(s.Digs);
            for (int r = 0; r < recs; r++) Advance(b, qb);
            rows.Add(("recompute_driver(timing,stability,settled,wake,fuzz)", b, true));
            var c = Make(s, "wet"); var sp = new Spreader(c, Recompute, DtCap);
            for (int t = 0; t < recs * Spreader.Ticks; t++) sp.Tick();
            rows.Add(("spread_driver(Spreader.Tick)", c, true));
            int edits = s.Digs.Count(d => d.t < T), fails = 0;
            foreach (var (name, m, wantZero) in rows)
            {
                bool ok = wantZero ? m.TailSteps <= edits : m.TailSteps > edits;
                if (!ok) fails++;
                Console.WriteLine((ok ? "ok   " : "FAIL ") + name + " tail_steps=" + m.TailSteps + " first_tail_t=" + m.FirstTailT.ToString("F3", Inv) + " steps=" + m.Steps + " t=" + m.T.ToString("F3", Inv));
            }
            Console.WriteLine((fails == 0 ? "EVENSTEPS PASS " + s.Name : "EVENSTEPS FAIL " + s.Name + " rows=" + fails) + " (edits in window: " + edits + ")");
            return fails == 0 ? 0 : 1;
        }

        // ------------------------------------------------------------------------------------------------ timing
        /// <summary>Advance one recompute (0.2 water-s) in CFL substeps, applying digs due. Returns substeps.</summary>
        static int Advance(LocalInertialDL m, Queue<(double t, bool[] mask, double[] depth)> digs)
        {
            while (digs != null && digs.Count > 0 && m.T >= digs.Peek().t - 1e-9) { var d = digs.Dequeue(); m.Dig(d.mask, d.depth); }
            return Spreader.AdvanceWhole(m, Recompute, DtCap);   // the stepping rule (LocalInertialDL.EvenDt)
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
                        var warm = Make(s, mode);
                        var wq = new Queue<(double, bool[], double[])>(s.Digs);
                        for (int k = 0; k < 25; k++) Advance(warm, wq);

                        var m = Make(s, mode);
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
        static int Stability(string dir, string outPath, string only, string mode)
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
                    var m = Make(s, mode);
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
                    string row = "{\"run\":\"" + file + "\",\"mode\":\"" + mode + "\",\"water_s\":" + F(seconds) + ",\"recomputes\":" + nrec + ",\"substeps\":" + sub
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
                    if ((k & 5) == 5) m.Sleep = true;          // a quarter: wet-only with settled water sleeping
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
                                    if (op == 0) m.RaiseBed(c, -(0.1 + 3 * rng.NextDouble()));       // dig
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

        // ---------------------------------------------------------------------------------------------- spread
        /// <summary>Each recompute spread over its 12 ticks: per-tick cost (mean, p95, max) against the unspread recompute,
        /// and the water checked bitwise against a twin run whole, at every recompute's end.</summary>
        static int SpreadTiming(string dir, string outPath, string[] modes, string only)
        {
            var onlySet = new HashSet<string>(only.Split(new[] { ',' }, StringSplitOptions.RemoveEmptyEntries));
            bool header = !File.Exists(outPath);
            int fails = 0;
            using (var w = new StreamWriter(outPath, append: true))
            {
                if (header) w.WriteLine("runtime,scene,n,mode,recomputes,whole_p95_ms_recompute,whole_max_ms_recompute,spread_mean_ms_tick,spread_p95_ms_tick,spread_max_ms_tick,identical,books_dl");
                foreach (var path in Directory.GetFiles(dir, "*.lis").OrderBy(p => p))
                {
                    var s = SceneFile.Read(path);
                    string scene = s.Name.Substring(0, s.Name.LastIndexOf('_'));
                    if (onlySet.Count > 0 && !onlySet.Contains(scene)) continue;
                    foreach (var mode in modes)
                    {
                        var warm = Make(s, mode); var wsp = new Spreader(warm, Recompute, DtCap);
                        for (int k = 0; k < 25 * Spreader.Ticks; k++) wsp.Tick();
                        var whole = Make(s, mode); var spread = Make(s, mode);
                        var sp = new Spreader(spread, Recompute, DtCap);
                        var q1 = new Queue<(double t, bool[] mask, double[] depth)>(s.Digs);
                        var q2 = new Queue<(double t, bool[] mask, double[] depth)>(s.Digs);
                        int nrec = (int)Math.Round(s.Duration / Recompute);
                        var tw = new double[nrec]; var tt = new double[nrec * Spreader.Ticks];
                        bool same = true;
                        for (int r = 0; r < nrec; r++)
                        {
                            long t0 = Stopwatch.GetTimestamp();
                            Advance(whole, q1);
                            tw[r] = (Stopwatch.GetTimestamp() - t0) * 1000.0 / Stopwatch.Frequency;
                            for (int k = 0; k < Spreader.Ticks; k++)
                            {
                                if (sp.AtRecomputeStart)
                                    while (q2.Count > 0 && spread.T >= q2.Peek().t - 1e-9) { var d = q2.Dequeue(); spread.Dig(d.mask, d.depth); }
                                t0 = Stopwatch.GetTimestamp();
                                sp.Tick();
                                tt[r * Spreader.Ticks + k] = (Stopwatch.GetTimestamp() - t0) * 1000.0 / Stopwatch.Frequency;
                            }
                            if (same && !whole.V.SequenceEqual(spread.V)) same = false;
                        }
                        var sw = tw.OrderBy(x => x).ToArray(); var st = tt.OrderBy(x => x).ToArray();
                        string row = string.Join(",", RuntimeLabel().Replace(',', ' '), scene, s.Nx.ToString(Inv), mode, nrec.ToString(Inv),
                            P95(sw).ToString("F4", Inv), sw[sw.Length - 1].ToString("F4", Inv),
                            tt.Average().ToString("F4", Inv), P95(st).ToString("F4", Inv), st[st.Length - 1].ToString("F4", Inv),
                            same ? "true" : "false", spread.BooksError.ToString(Inv));
                        w.WriteLine(row); w.Flush(); Console.WriteLine(row);
                        if (!same || spread.BooksError != 0) fails++;
                    }
                }
            }
            Console.WriteLine(fails == 0 ? "SPREAD PASS" : "SPREAD FAIL rows=" + fails);
            return fails == 0 ? 0 : 1;
        }

        static double P95(double[] sorted) => sorted[(int)Math.Min(sorted.Length - 1, Math.Ceiling(0.95 * sorted.Length) - 1)];

        // --------------------------------------------------------------------------------------------- settled
        /// <summary>Settled water sleeping: run until every wet block sleeps (cap 100 water-s), then time 300 recomputes
        /// of settled water (per-tick = recompute / 12), and check it stayed put.</summary>
        static int Settled(string dir, string outPath, string only)
        {
            var onlySet = new HashSet<string>(only.Split(new[] { ',' }, StringSplitOptions.RemoveEmptyEntries));
            bool header = !File.Exists(outPath);
            using (var w = new StreamWriter(outPath, append: true))
            {
                if (header) w.WriteLine("runtime,scene,n,recomputes_to_sleep,wet_blocks,sleeping_blocks,mean_ms_tick,p95_ms_tick,max_ms_recompute,surface_move_m,books_dl");
                foreach (var path in Directory.GetFiles(dir, "*.lis").OrderBy(p => p))
                {
                    var s = SceneFile.Read(path);
                    string scene = s.Name.Substring(0, s.Name.LastIndexOf('_'));
                    if (onlySet.Count > 0 && !onlySet.Contains(scene)) continue;
                    var m = Make(s, "sleep");
                    int toSleep = -1, cap = (int)(100 / Recompute);
                    for (int r = 0; r < cap; r++)
                    {
                        Advance(m, null);
                        if (m.SleepingBlocks == WetBlocks(m)) { toSleep = r + 1; break; }
                    }
                    var eta0 = Eta(m);
                    var ms = new double[300];
                    for (int r = 0; r < ms.Length; r++)
                    {
                        long t0 = Stopwatch.GetTimestamp();
                        Advance(m, null);
                        ms[r] = (Stopwatch.GetTimestamp() - t0) * 1000.0 / Stopwatch.Frequency;
                    }
                    double mv = 0; var e = Eta(m);
                    for (int c = 0; c < e.Length; c++) if (m.V[c] > 0) mv = Math.Max(mv, Math.Abs(e[c] - eta0[c]));
                    var so = ms.OrderBy(x => x).ToArray();
                    string row = string.Join(",", RuntimeLabel().Replace(',', ' '), scene, s.Nx.ToString(Inv), toSleep.ToString(Inv),
                        WetBlocks(m).ToString(Inv), m.SleepingBlocks.ToString(Inv), (ms.Average() / 12).ToString("F5", Inv),
                        (P95(so) / 12).ToString("F5", Inv), so[so.Length - 1].ToString("F4", Inv), mv.ToString("E2", Inv), m.BooksError.ToString(Inv));
                    w.WriteLine(row); w.Flush(); Console.WriteLine(row);
                }
            }
            return 0;
        }

        static int WetBlocks(LocalInertialDL m)
        {
            int nbx = (m.Nx + LocalInertialDL.B - 1) / LocalInertialDL.B, nby = (m.Ny + LocalInertialDL.B - 1) / LocalInertialDL.B, n = 0;
            for (int bj = 0; bj < nby; bj++)
                for (int bi = 0; bi < nbx; bi++)
                {
                    bool wet = false;
                    for (int j = bj * LocalInertialDL.B; j < Math.Min((bj + 1) * LocalInertialDL.B, m.Ny) && !wet; j++)
                        for (int i = bi * LocalInertialDL.B; i < Math.Min((bi + 1) * LocalInertialDL.B, m.Nx); i++)
                            if (m.V[j * m.Nx + i] > 0) { wet = true; break; }
                    if (wet) n++;
                }
            return n;
        }

        // ------------------------------------------------------------------------------------------------ wake
        /// <summary>Wake check: let the scene's water fall asleep (sleep mode), with a never-sleeping twin (wet mode) run
        /// alongside; then dig a 24x24 hole 2 m deep under the water's deepest square and run 20 water-s more. The sleeping
        /// model must match the twin within the fidelity tolerances (L1 2%, IoU 0.98) and the hole's mean depth within 2%.
        /// LITIMING_MUTATE=nowake (edits do not wake) must turn it red.</summary>
        static int Wake(string file, string outPath)
        {
            var s = SceneFile.Read(file);
            var slp = Make(s, "sleep"); var twin = Make(s, "wet");
            int r0 = 0;
            for (; r0 < (int)(100 / Recompute); r0++)
            {
                Advance(slp, null); Advance(twin, null);
                if (slp.SleepingBlocks == WetBlocks(slp)) break;
            }
            int asleep = slp.SleepingBlocks, deep = 0;
            for (int c = 0; c < slp.V.Length; c++) if (slp.V[c] > slp.V[deep]) deep = c;
            var mask = new bool[slp.V.Length]; var depth = new double[slp.V.Length];
            int dj = deep / s.Nx, di = deep % s.Nx;
            for (int j = Math.Max(0, dj - 12); j < Math.Min(s.Ny, dj + 12); j++)
                for (int i = Math.Max(0, di - 12); i < Math.Min(s.Nx, di + 12); i++) { mask[j * s.Nx + i] = true; depth[j * s.Nx + i] = 2.0; }
            slp.Dig(mask, depth); twin.Dig(mask, depth);
            for (int r = 0; r < (int)(20 / Recompute); r++) { Advance(slp, null); Advance(twin, null); }
            var (l1, iou) = Compare(slp.DepthArray(), twin.DepthArray());
            var hs = slp.DepthArray(); var ht = twin.DepthArray();
            double a = 0, b = 0; int n = 0;
            for (int c = 0; c < hs.Length; c++) if (mask[c]) { a += hs[c]; b += ht[c]; n++; }
            // the hole must have filled: its surface back near the lake's (a stuck hole keeps the 2 m drop)
            double holeRel = Math.Abs(a - b) / b;
            bool ok = l1 <= 0.02 && iou >= 0.98 && holeRel <= 0.02 && slp.BooksError == 0 && asleep > 0;
            // and the disturbed water settles again: recomputes until every wet block sleeps once more (cap 600 water-s)
            int resleep = -1;
            double resleepCap = double.Parse(Environment.GetEnvironmentVariable("LITIMING_RESLEEP_S") ?? "600", Inv);
            bool trace = Environment.GetEnvironmentVariable("LITIMING_RESLEEP_S") != null;
            for (int r = 0; r < (int)(resleepCap / Recompute); r++)
            {
                Advance(slp, null);
                if (trace && r % 500 == 0) Console.WriteLine("resleep t=" + (r * Recompute).ToString("F0", Inv) + " sleeping " + slp.SleepingBlocks + "/" + WetBlocks(slp) + " active " + slp.ActiveBlocks);
                if (slp.SleepingBlocks == WetBlocks(slp)) { resleep = r + 1; break; }
            }
            string row = "{\"file\":\"" + Path.GetFileName(file) + "\",\"recomputes_to_sleep\":" + (r0 + 1) + ",\"asleep_blocks_at_dig\":" + asleep
                + ",\"l1rel\":" + F(l1) + ",\"iou\":" + F(iou) + ",\"hole_mean_sleep_m\":" + F(a / n) + ",\"hole_mean_twin_m\":" + F(b / n)
                + ",\"hole_rel_err\":" + F(holeRel) + ",\"recomputes_to_resleep\":" + resleep + ",\"sleeping_after\":" + slp.SleepingBlocks
                + ",\"wet_blocks_after\":" + WetBlocks(slp) + ",\"books_dl\":" + slp.BooksError + ",\"pass\":" + (ok ? "true" : "false") + "}";
            File.AppendAllText(outPath, row + Environment.NewLine);
            Console.WriteLine(row);
            Console.WriteLine(ok ? "WAKE PASS" : "WAKE FAIL");
            return ok ? 0 : 1;
        }
    }
}
