using System;
using System.Collections.Generic;
using System.IO;

namespace LiTiming
{
    /// <summary>A scene exported by <c>mandrakemodels.csharp_export</c> (format "LIS2", little-endian; "LIS1" lacks FrMax and
    /// reads as uncapped).</summary>
    public sealed class SceneFile
    {
        public string Name;
        public int Ny, Nx;
        public double Dx, Manning, Theta, Alpha, Duration;
        /// <summary>The reference solver's face Froude cap; 0 = none (plain local_inertial).</summary>
        public double FrMax;
        public double[] Z, H;
        public bool[] Wall;
        public double Evaporation, Rain;
        public List<(bool[] mask, double q)> Sources = new List<(bool[], double)>();
        public List<(bool[] mask, double rate)> Sinks = new List<(bool[], double)>();
        public List<(bool[] mask, int kind, double a, double b)> Levels = new List<(bool[], int, double, double)>();
        public bool HasGroundwater;
        public double GwTable, GwK;
        public List<(double t, bool[] mask, double[] depth)> Digs = new List<(double, bool[], double[])>();
        public List<(double t, int steps, double verr, double[] h)> Checkpoints = new List<(double, int, double, double[])>();

        public static SceneFile Read(string path)
        {
            var s = new SceneFile { Name = Path.GetFileNameWithoutExtension(path) };
            using (var r = new BinaryReader(File.OpenRead(path)))
            {
                var magic = new string(r.ReadChars(4));
                if (magic != "LIS1" && magic != "LIS2") throw new InvalidDataException(path + ": not a LIS1/LIS2 file");
                s.Ny = r.ReadInt32(); s.Nx = r.ReadInt32();
                s.Dx = r.ReadDouble(); s.Manning = r.ReadDouble(); s.Theta = r.ReadDouble(); s.Alpha = r.ReadDouble();
                if (magic == "LIS2") s.FrMax = r.ReadDouble();
                s.Duration = r.ReadDouble();
                int n = s.Ny * s.Nx;
                s.Z = Doubles(r, n); s.H = Doubles(r, n); s.Wall = Mask(r, n);
                s.Evaporation = r.ReadDouble(); s.Rain = r.ReadDouble();
                for (int k = r.ReadInt32(); k > 0; k--) { var m = Mask(r, n); s.Sources.Add((m, r.ReadDouble())); }
                for (int k = r.ReadInt32(); k > 0; k--) { var m = Mask(r, n); s.Sinks.Add((m, r.ReadDouble())); }
                for (int k = r.ReadInt32(); k > 0; k--)
                {
                    var m = Mask(r, n); int kind = r.ReadInt32(); double a = r.ReadDouble(), b = r.ReadDouble();
                    s.Levels.Add((m, kind, a, b));
                }
                s.HasGroundwater = r.ReadInt32() != 0;
                if (s.HasGroundwater) { s.GwTable = r.ReadDouble(); s.GwK = r.ReadDouble(); }
                for (int k = r.ReadInt32(); k > 0; k--)
                {
                    double t = r.ReadDouble(); var m = Mask(r, n); s.Digs.Add((t, m, Doubles(r, n)));
                }
                for (int k = r.ReadInt32(); k > 0; k--)
                {
                    double t = r.ReadDouble(); int steps = r.ReadInt32(); double verr = r.ReadDouble();
                    s.Checkpoints.Add((t, steps, verr, Doubles(r, n)));
                }
            }
            return s;
        }

        static double[] Doubles(BinaryReader r, int n)
        {
            var bytes = r.ReadBytes(n * 8);
            var a = new double[n];
            Buffer.BlockCopy(bytes, 0, a, 0, bytes.Length);
            return a;
        }

        static bool[] Mask(BinaryReader r, int n)
        {
            var bytes = r.ReadBytes(n);
            var m = new bool[n];
            for (int i = 0; i < n; i++) m[i] = bytes[i] != 0;
            return m;
        }
    }
}
