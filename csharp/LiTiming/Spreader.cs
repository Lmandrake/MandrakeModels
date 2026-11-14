using System;

namespace LiTiming
{
    /// <summary>
    /// Spreads each 12-tick recompute's work over its 12 ticks (stage C finding: unspread, the whole recompute lands in one
    /// tick, 17 to 20 ms on the 325 stress load). A recompute is k equal substeps (k = ceil(0.2 s / CFL dt), fixed at the
    /// recompute's start: a ragged split makes local-inertial slosh); each substep is four passes over the active blocks
    /// (faces, limit, move, then the step's end, counted as one more pass). Each tick does its share of the work left,
    /// counted in block-passes, so the last tick of a recompute always finishes it.
    ///
    /// Water is bitwise the same as running the recompute in one go (the harness checks it): pass 1 reads only the
    /// previous state, and a pass is never interleaved with another substep. Between ticks the water is mid-step
    /// (some faces moved, some not); every move is whole and booked, so the books are exact at every tick. Edits are taken
    /// only at a recompute's start (the game would queue a player's dig until then).
    /// </summary>
    public sealed class Spreader
    {
        public const int Ticks = 12;
        readonly LocalInertialDL m;
        readonly double recompute, dtCap;
        int tick = Ticks;            // ticks done in the current recompute (Ticks = none running)
        int k, sub, phase, cursor;   // substeps this recompute, substep index, pass 0..3, next active block in the pass
        double dt;

        public Spreader(LocalInertialDL m, double recompute = 0.2, double dtCap = 5.0)
        {
            this.m = m; this.recompute = recompute; this.dtCap = dtCap;
        }

        /// <summary>Substeps in the recompute running now (or the last one).</summary>
        public int Substeps => k;
        public bool AtRecomputeStart => tick == Ticks;

        static int SubstepsFor(LocalInertialDL m, double recompute, double dtCap)
        {
            int k = (int)Math.Ceiling(recompute / Math.Min(m.StableDt(), dtCap) - 1e-9);
            return k < 1 ? 1 : k;
        }

        /// <summary>The whole recompute in one go (the unspread reference): k equal substeps.</summary>
        public static int AdvanceWhole(LocalInertialDL m, double recompute = 0.2, double dtCap = 5.0)
        {
            int k = SubstepsFor(m, recompute, dtCap);
            double dt = recompute / k;
            for (int s = 0; s < k; s++) m.Step(dt);
            return k;
        }

        /// <summary>One game tick of water work. Returns true on the tick that finishes a recompute.</summary>
        public bool Tick()
        {
            if (tick == Ticks)
            {
                tick = 0;
                k = SubstepsFor(m, recompute, dtCap);
                dt = recompute / k;
                sub = 0; phase = 0; cursor = 0;
                m.BeginStep(dt);
            }
            int ticksLeft = Ticks - tick;
            tick++;
            if (sub >= k) return tick == Ticks;          // finished early (cannot happen: the last tick takes the rest)
            // work left, in block-passes; later substeps are estimated at today's active count
            long n = m.ActiveBlocks;
            long left = (n * (4 - phase) - cursor) + (long)(k - sub - 1) * 4 * Math.Max(n, 1);
            long budget = ticksLeft <= 1 ? long.MaxValue : (left + ticksLeft - 1) / ticksLeft;
            while (budget > 0 && sub < k)
            {
                n = m.ActiveBlocks;
                if (phase == 3)
                {
                    m.EndStep();
                    budget -= Math.Max(n, 1);
                    phase = 0; cursor = 0; sub++;
                    if (sub < k) m.BeginStep(dt);
                    continue;
                }
                int take = (int)Math.Min(budget, n - cursor);
                int a1 = cursor + take;
                if (phase == 0) m.FacesRange(cursor, a1);
                else if (phase == 1) m.LimitRange(cursor, a1);
                else m.MoveRange(cursor, a1);
                budget -= Math.Max(take, 1);
                cursor = a1;
                if (cursor >= n) { phase++; cursor = 0; }
            }
            return tick == Ticks;
        }
    }
}
