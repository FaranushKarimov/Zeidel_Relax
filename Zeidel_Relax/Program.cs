using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
namespace Zeidel_Relax
{
    class Program
    {
        static int n = 20;
        static double h = 1.0 / n;
        static double eps = Math.Pow(h, 5);
        static double[] x = new double[n + 1];
        static double[] ai = new double[n + 1];
        static double[] gi = new double[n + 1];
        static double[] fi = new double[n + 1];
        static double[] alpha = new double[n + 1];
        static double[] beta = new double[n + 1];
        static double[] a_big = new double[n + 1];
        static double[] b_big = new double[n + 1];
        static double[] c_big = new double[n + 1];
        static double[,] A = new double[3, n + 1];
        static double[] y = new double[n + 1];
        static double[] y_zeidel = new double[n + 1];
        static double[] y_relax = new double[n + 1];
        static double[] y_grad = new double[n + 1];
        static int k;
        static int best_k = 99999;
        static double best_omega;
        public static double Px(double x)
        {
            return 1 + Math.Pow(x, 3);
        }
        public static double Gx(double x)
        {
            return x + 1;
        }
        public static double Ux(double x)
        {
            return Math.Pow(x,3) * (1 - x);
        }
        public static double Fx(double x)
        {
            return 23 * Math.Pow(x, 5) - 15 * Math.Pow(x, 4) + Math.Pow(x, 3) + 12 * Math.Pow(x, 2) - 6 * x;
        }

        static void Main(string[] args)
        {
            alpha[1] = 0;
            beta[1] = 0;
            y[0] = 0;
            y[n] = 0;
            // вычисление a[i],g[i],f[i]
            for (int i = 0; i <= n; i++)
            {
                ai[i] = Px(i * h);
                gi[i] = Gx(i * h);
                fi[i] = Fx(i * h) * Math.Pow(h, 2);
            }
            for (int i = 0; i < n; i++)
                x[i] = i * h;
            // вычисление A,B,C
            for (int i = 1; i < n; i++)
            {
                a_big[i] = -ai[i];
                b_big[i] = -(ai[i] + ai[i + 1] + h * h * gi[i]);
                c_big[i] = -ai[i + 1];
            }
            for (int i = 0; i < n; i++)
                A[0, i] = a_big[i];
            for (int i = 0; i < n; i++)
                A[1, i] = b_big[i];
            for (int i = 0; i < n; i++)
                A[2, i] = c_big[i];
            // вычисление alpha, beta
            for (int i = 1; i < n; i++)
            {
                alpha[i + 1] = c_big[i] / (b_big[i] - a_big[i] * alpha[i]);
                beta[i + 1] = (a_big[i] * beta[i] - fi[i]) / (b_big[i] - a_big[i] * alpha[i]);
            }
            // вычисление y[i]
            for (int i = n - 1; i > 0; i--)
            {
                y[i] = alpha[i + 1] * y[i + 1] + beta[i + 1];
            }
            y_zeidel = Zeidel(ai, fi, gi);
            //вывод метода Зейделя
            using (StreamWriter sw = new StreamWriter("C:/Users/User/source/repos/Zeidel_Relax/Zeidel.txt"))
            {
                for (int i = 0; i <= n; i++)
                    sw.WriteLine("{0}\t{1}\t{2}\t{3:0.000000000000}", i * h, y[i], y_zeidel[i], Math.Abs(y[i] - y_zeidel[i]));
                sw.WriteLine("Количество итераций: {0}", k);
            }
            y_relax = Relax(ai, fi, gi);
            //вывод метода релаксации
            using (StreamWriter sw = new StreamWriter("C:/Users/User/source/repos/Zeidel_Relax/Relax.txt"))
            {
                sw.WriteLine("Вычисления для омега={0}", best_omega);
                for (int i = 0; i <= n; i++)
                    sw.WriteLine("{0}\t{1}\t{2}\t{3:0.000000000}", i * h, y[i], y_relax[i], Math.Abs(y[i] - y_relax[i]));
            }
            y_grad = Gradient(ai, fi, gi, n, h, eps);
            //вывод метода градиентного спуска
            using (StreamWriter sw = new StreamWriter("C:/Users/User/source/repos/Zeidel_Relax/Gradient.txt"))
            {
                for (int i = 0; i <= n; i++)
                    sw.WriteLine("{0}\t{1}\t{2}\t{3:0.000000000}", i * h, y[i], y_grad[i], Math.Abs(y[i] - y_grad[i]));
            }
        }
        public static double[] Zeidel(double[] a, double[] f, double[] g)
        {
            double[] y = new double[n + 1];
            double r;
            k = 0;
            for (int i = 0; i < n; i++)
                y[i] = 0;
            do
            {
                r = -1;
                for (int i = 1; i < n; i++)
                    y[i] = (f[i] + a[i + 1] * y[i + 1] + a[i] * y[i - 1]) / (a[i] + a[i + 1] + g[i] * h * h);
                k++;
                for (int i = 1; i < n; i++)
                    if (Math.Abs(-a[i + 1] * y[i + 1] + (a[i + 1] + a[i] + g[i] * h * h) * y[i] - a[i] * y[i - 1] - f[i]) > r)
                        r = Math.Abs(-a[i + 1] * y[i + 1] + (a[i + 1] + a[i] + g[i] * h * h) * y[i] - a[i] * y[i - 1] - f[i]);
            } while (r > eps);
            return y;
        }
        public static double[] Relax(double[] a, double[] f, double[] g)
        {
            double[] y = new double[n + 1];
            double[] y_pre = new double[n + 1];
            double r = 0;
            for (double omega = 1.05; omega < 2; omega += 0.05)
            {
                for (int i = 0; i <= n; i++)
                    y[i] = 0;
                k = 0;
                do
                {
                    for (int i = 0; i <= n; i++)
                        y_pre[i] = y[i];
                    r = -1;
                    for (int i = 1; i < n; i++)
                        y[i] = (f[i] + a[i + 1] * y_pre[i + 1] + a[i] * y[i - 1]) / (a[i + 1] + a[i] + g[i] * h * h) * omega + (1 - omega) * y_pre[i];
                    for (int i = 1; i < n; i++)
                        if (Math.Abs(-a[i + 1] * y[i + 1] + (a[i + 1] + a[i] + g[i] * h * h) * y[i] - a[i] * y[i - 1] - f[i]) > r)
                            r = Math.Abs(-a[i + 1] * y[i + 1] + (a[i + 1] + a[i] + g[i] * h * h) * y[i] - a[i] * y[i - 1] - f[i]);
                    k++;
                } while (r > eps);
                if (k < best_k)
                {
                    best_omega = omega;
                    best_k = k;
                }
                using (StreamWriter sw = new StreamWriter("C:/Users/User/source/repos/Zeidel_Relax/omega.txt", true))
                {
                    sw.WriteLine("{0}\t{1}\n", omega, k);
                }
            }
            for (int i = 0; i <= n; i++)
                y[i] = 0;
            k = 0;
            do
            {
                for (int i = 0; i <= n; i++)
                    y_pre[i] = y[i];
                r = -1;
                for (int i = 1; i < n; i++)
                    y[i] = (f[i] + a[i + 1] * y_pre[i + 1] + a[i] * y[i - 1]) / (a[i + 1] + a[i] + g[i] * h * h) * best_omega + (1 - best_omega) * y_pre[i];
                k++;
                for (int i = 1; i < n; i++)
                    if (Math.Abs(-a[i + 1] * y[i + 1] + (a[i + 1] + a[i] + g[i] * h * h) * y[i] - a[i] * y[i - 1] - f[i]) > r)
                        r = Math.Abs(-a[i + 1] * y[i + 1] + (a[i + 1] + a[i] + g[i] * h * h) * y[i] - a[i] * y[i - 1] - f[i]);
            } while (r > eps);
            return y;
        }
        public static double[] Gradient(double[] a, double[] f, double[] g, int n, double h, double eps)
        {
            double[] y = new double[n + 1];
            double[] y_pre = new double[n + 1];
            double r = 0;
            for (int i = 0; i <= n; i++)
                y[i] = 0;
            int k = 0;
            double[] rr = new double[n + 1];
            for (int i = 0; i <= n; i++)
                rr[i] = 0;
            double sk_1;
            double sk_2;
            double res;
            do
            {
                for (int i = 0; i <= n; i++)
                    y_pre[i] = y[i];
                r = -1;
                sk_1 = 0;
                sk_2 = 0;
                res = 0;
                for (int i = 1; i < n; i++)
                    sk_1 = sk_1 + rr[i] * rr[i];
                for (int i = 1; i < n; i++)
                    sk_2 = sk_2 + rr[i] * (-a[i + 1] * rr[i + 1] + (a[i + 1] + a[i] + g[i] * h * h) * rr[i] - a[i] * rr[i - 1]);
                if (sk_2 != 0)
                    res = sk_1 / sk_2;
                for (int i = 1; i < n; i++)
                {
                    y[i] = y_pre[i] - res * rr[i];
                }
                k++;
                for (int i = 1; i < n; i++)
                    rr[i] = -a[i + 1] * y[i + 1] + (a[i + 1] + a[i] + g[i] * h * h) * y[i] - a[i] * y[i - 1] - f[i];
                for (int i = 1; i < n; i++)
                    if (Math.Abs(-a[i + 1] * y[i + 1] + (a[i + 1] + a[i] + g[i] * h * h) * y[i] - a[i] * y[i - 1] - f[i]) > r)
                        r = Math.Abs(-a[i + 1] * y[i + 1] + (a[i + 1] + a[i] + g[i] * h * h) * y[i] - a[i] * y[i - 1] - f[i]);
            } while (r > eps);

            using (StreamWriter sw = new StreamWriter("C:/Users/User/source/repos/Zeidel_Relax/number_of_iter.txt", true))
            {
                sw.WriteLine("{0}\n", k);
            }
            return y;
        }
    }
}
    
