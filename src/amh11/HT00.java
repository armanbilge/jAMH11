/**
 * HT00.java
 * 
 * AMH11: Java implementation of the matrix exponential method
 *     described by Al-Mohy and Higham (2011)
 * 
 * Copyright (C) 2014 Arman D. Bilge <armanbilge@gmail.com>
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package amh11;

import java.util.Arrays;
import java.util.Comparator;
import java.util.Random;

import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;

/**
 * @author Arman D. Bilge <armanbilge@gmail.com>
 *
 */
public final class HT00 {

    private HT00() {}
    
    public static interface HT00Random {
        public boolean nextBoolean();
        public double nextDouble();
    }
    
    private static HT00Random random = new HT00Random() {
        private final Random random = new Random();
        public boolean nextBoolean() {
            return random.nextBoolean();
        }
        public double nextDouble() {
            return random.nextDouble();
        }
    };
    
    public static final void setRandom(HT00Random random) {
        HT00.random = random;
    }
        
    public static interface MatrixFunction {
        public Matrix apply(Matrix X, boolean transpose);
        public int getDimensions();
        public boolean isReal();
    }
    
    public static final double[] normest1(MatrixFunction f, int t) {
        
        int n = f.getDimensions();
        boolean isReal = f.isReal();
        
        t = Math.abs(t);
        
//        int rpt_S = 0;
//        int rpt_e = 0;
        
        if (t == n || n <= 4) {
            Matrix X = DoubleFactory2D.dense.identity(n);
            Matrix Y = f.apply(X, false);
            
            Vector sums = absColumnSums(Y);
            double est = Algebra.DEFAULT.normInfinity(sums);
            return new double[]{est, 1};
        }
        
        Matrix X = randomSigns(n, t);
        for (int i = 0; i < n; ++i) X.setQuick(i, 0, 1.0);
        unduplicate(X, null);
        X.assign(Functions.div(n));
        
        int itMax = 5;
        int it = 0;
        int nmv = 0;
        
        Vector ind = DoubleFactory1D.dense.make(t);
        Vector vals = DoubleFactory1D.dense.make(t);
        Vector Zvals = DoubleFactory1D.dense.make(n);
        Matrix S = DoubleFactory2D.dense.make(n, t);
        Matrix S_old;
        Vector indHist = null;
        
        double est;
        double estOld = 0.0;
        
        while (true) {
            
            ++it;
            Matrix Y = f.apply(X, false);
            ++nmv;
            
            vals = absColumnSums(Y);
            final double[] valsArray = vals.toArray();
            Integer[] m = new Integer[vals.size()];
            for (int i = 0; i < m.length; ++i) m[i] = i;
            Arrays.sort(m, new Comparator<Integer>() {
                public int compare(Integer i, Integer j) {
                    return Double.compare(valsArray[i], valsArray[j]);
                }
            });
            int[] m2 = new int[t];
            for (int i = 0; i < t; ++i) {
                m2[i] = m[t-1-i];
            }
            vals = vals.viewSelection(m2);
            Vector vals_ind = ind.viewSelection(m2);
            est = vals.getQuick(0);
            
            int est_j = Integer.MIN_VALUE;
            if (est > estOld || (it == 2)) {
                est_j = (int) vals_ind.getQuick(0);
            }
            
            if (it >= 2 && est <= estOld) {
                est = estOld;
                break;
            }
            estOld = est;
            
            if (it > itMax) {
                it = itMax;
                break;
            }
            
            S_old = S;
            S = mySign(Y);
            
            if (isReal) {
                Matrix SS = Algebra.DEFAULT
                        .mult(Algebra.DEFAULT.transpose(S_old), S);
                double np = columnsMax(SS).zSum();
                if (np == t) break;
                /* double r = */ unduplicate(S, S_old);
//                rpt_S += r;
            }
            
            Matrix Z = f.apply(S, true);
            ++nmv;
            for (int i = 0; i < n; ++i)
                Zvals.setQuick(i, Algebra.DEFAULT.normInfinity(Z.viewRow(i)));
            
            if (it >= 2 && Algebra.DEFAULT.normInfinity(Zvals)
                        == Zvals.getQuick(est_j)) break;
            
            m = new Integer[Zvals.size()];
            final double[] zvalsArray = Zvals.toArray();
            for (int i = 0; i < m.length; ++i) m[i] = i;
            Arrays.sort(m, new Comparator<Integer>() {
                public int compare(Integer i, Integer j) {
                    return Double.compare(zvalsArray[i], zvalsArray[j]);
                }
            });
            m2 = new int[n];
            for (int i = 0; i < n; ++i) m2[i] = m[n-1-i];
            
            int imax = t;
            int[] m2Tot = Arrays.copyOfRange(m2, 0, t);
            if (it == 1) {
                ind = DoubleFactory1D.dense.make(t);
                for (int i = 0; i < t; ++i) ind.setQuick(i, m2Tot[i]);
                indHist = ind;
            } else {
                int rep = memberCount(m2Tot, indHist);
//                rpt_e += rep;
                if (rep == t) break;
                int j = 0;
                for (int i = 0; i < t; ++i) {
                    if (j >= n) {
                        imax = i-1;
                        break;
                    }
                    while (contains(indHist, m2[j])) {
                        ++j;
                        if (j >= n) {
                            imax = i - 1;
                            break;
                        }
                    }
                    if (j >= n) break;
                    ind.setQuick(i, m2[j]);
                    ++j;
                }
                Vector temp = DoubleFactory1D.dense.make(t + imax);
                for (int i = 0; i < t; ++i)
                    temp.setQuick(i, indHist.getQuick(i));
                for (int i = 0; i < imax; ++i)
                    temp.setQuick(t+i, ind.getQuick(i));
                indHist = temp;
            }
            X = DoubleFactory2D.dense.make(n, t);
            for (int j = 0; j < imax; ++j)
                X.setQuick((int) ind.getQuick(j), j, 1);
        }
        
        int iter = nmv;
        return new double[]{est, iter};
                
    }
    
    
    private static final Vector absColumnSums(Matrix X) {
        Vector sums = DoubleFactory1D.dense.make(X.columns());
        for (int i = 0; i < X.columns(); ++i) {
            double sum = Algebra.DEFAULT.norm1(X.viewColumn(i));
            sums.setQuick(i, sum);
        }
        return sums;
    }
    
    private static final Matrix randomSigns(int r, int c) {
        Matrix signs = DoubleFactory2D.dense.make(r, c);
        for (int i = 0; i < r; ++i) {
            for (int j = 0; j < c; ++j) {
                double sign = random.nextBoolean() ? 1.0 : -1.0;
                signs.setQuick(i, j, sign);
            }
        }
        return signs;
    }
    
    private static final double unduplicate(Matrix S,
            Matrix S_old) {
        
        int n = S.rows();
        int t = S.columns();
        double r = 0.0;
        if (t == 1) return r;
        
        Matrix W;
        int jStart;
        int lastCol;
        if (S_old == null) {
            W = DoubleFactory2D.dense.make(n, t);
            for (int i = 0; i < n; ++i)
                W.setQuick(i, 0, S.getQuick(i, 0));
            jStart = 1;
            lastCol = 1;
        } else {
            W = DoubleFactory2D.dense.make(n, 2*t-1);
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < t; ++j) {
                    W.setQuick(i, j, S_old.getQuick(i, j));
                }
            }
            jStart = 0;
            lastCol = t;
        }
        
        for (int j = jStart; j < t; ++j) {
            int rpt = 0;
            while (Algebra.DEFAULT.norm2(Algebra.DEFAULT.mult(
                    Algebra.DEFAULT.transpose(S.viewPart(0, j, n, 1)),
                    W.viewPart(0, 0, n, lastCol)).assign(Functions.abs)) == n) {
                ++rpt;
                Matrix signs = randomSigns(n, 1);
                for (int i = 0; i < n; ++ i)
                    S.setQuick(i, j, signs.getQuick(i, 0));
                if (rpt > n / (double) t) break;                
            }
            
            r += Math.signum(rpt);
            if (j < t) {
                ++lastCol;
                for (int i = 0; i < n; ++i) {
                    W.setQuick(i, lastCol-1, S.getQuick(i, j));
                }
            }
        }
        
        return r;
    }
        
    private static final Matrix mySign(Matrix X) {
        Matrix signs =
                DoubleFactory2D.dense.make(X.rows(), X.columns());
        for (int i = 0; i < X.rows(); ++i) {
            for (int j = 0; j < X.columns(); ++j) {
                double sign = Math.signum(X.getQuick(i, j));
                if (sign == 0) sign = 1.0;
                signs.setQuick(i, j, sign);
            }
        }
        return signs;
    }
    
    private static final Vector columnsMax(Matrix X) {
        Vector maxes = DoubleFactory1D.dense.make(X.columns());
        for (int i = 0; i < maxes.size(); ++i)
            maxes.setQuick(i, Algebra.DEFAULT.normInfinity(X.viewColumn(i)));
        return maxes;
    }
    
    private static final int memberCount(int[] a, Vector b) {
        int count = 0;
        for (int i = 0; i < a.length; ++i) {
            double d = a[i];
            if (contains(b, d)) ++count;
        }
        return count;
    }
    
    private static final boolean contains(Vector v, double d) {
        for (int i = 0; i < v.size(); ++i) {
            if (v.getQuick(i) == d) return true;
        }
        return false;
    }
    
}
