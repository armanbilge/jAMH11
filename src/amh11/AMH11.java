/**
 * AMH11.java
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

import org.jblas.DoubleMatrix;

import amh11.HT00.DoubleMatrixFunction;

/**
 * @author Arman D. Bilge <armanbilge@gmail.com>
 *
 */
public final class AMH11 {
    
    private static final double TOL = Math.pow(2, -53);
    
    private AMH11() {}
    
    public static final DoubleMatrix expmv(double t, DoubleMatrix A,
            DoubleMatrix b) {
        return expmv(t, A, b, null, true, false, true);
    }
    
    public static final DoubleMatrix expmv(double t, DoubleMatrix A,
            DoubleMatrix b, DoubleMatrix M, boolean shift, boolean bal,
            boolean fullTerm) {
                
        A = A.dup();
        b = b.dup();
        
        if (bal) {
            throw new RuntimeException("Not implemented!");
        }
        
        int n = A.getRows();
        double mu = 0.0;
        if (shift) {
            mu = Utils.trace(A) / n;
            A.mini(DoubleMatrix.eye(n).muli(mu));
        }
        double tt;
        if (M == null) {
            tt = 1.0;
            M = selectTaylorDegree(A.mul(t), b, 55, 8, shift, bal, false);
        } else {
            tt = t;
        }
        
        double s = 1.0;
        int m, mMax;
        if (t == 0) {
            m = 0;
        } else {
            mMax = M.getRows();
            DoubleMatrix U = Utils.ascendingMatrix(mMax);
            DoubleMatrix C =
                    Utils.ceil(M.muli(Math.abs(tt))).transpose().mmuli(U);
            Utils.Zero2Inf(C);
            double cost;
            int[] min = new int[2];
            cost = getMin(C, min);
            m = min[C.getColumns() > 1 ? 1 : 0] + 1;
            if (cost == Double.POSITIVE_INFINITY)
                cost = 0;
            s = Math.max(cost/m, 1);
        }
        
        double eta = 1;
        if (shift) eta = Math.exp(t * mu / s);
        DoubleMatrix f = b.dup();
        for (int i = 0; i < s; ++i) {
            double c1 = b.normmax();
            for (int k = 1; k <= m; ++k) {
                b = A.mmul(b).mul(t/(s*k));
                f.addi(b);
                double c2 = b.normmax();
                if (!fullTerm) {
                    if (c1 + c2 <= TOL * f.normmax()) break;
                    c1 = c2;
                }
            }
            b = f.muli(eta);
        }
        return f;
    }
        
    private static final double getMin(DoubleMatrix M, int[] index) {
        double min = Double.POSITIVE_INFINITY;
        for (int i = 0; i < M.getRows(); ++i) {
            for (int j = 0; j < M.getColumns(); ++j) {
                double d = M.get(i, j);
                if (d < min) {
                    min = d;
                    index[0] = i;
                    index[1] = j;
                }
            }
        }
        return min;
    }
    
    public static final DoubleMatrix selectTaylorDegree(DoubleMatrix A,
            DoubleMatrix b, int mMax, int pMax, boolean shift,
            boolean bal, boolean forceEstm) {
        int n = A.getRows();
        if (bal) {
            throw new RuntimeException("Not implemented!");
        }
        double mu;
        if (shift) {
            mu = Utils.trace(A) / n;
            A.addi(DoubleMatrix.eye(n).muli(-mu));
        }
        double mv = 0.0;
        double normA = 0.0;
        if (!forceEstm)
            normA = Algebra.DEFAULT.norm1(A);
        double[] alpha;
        if (!forceEstm && normA < 4 * ThetaTaylor.THETA[mMax] * pMax
                * (pMax + 3) / (mMax * b.size())) {
            alpha = new double[pMax - 1];
            Arrays.fill(alpha, normA);
        } else {
            double[] eta = new double[pMax];
            alpha = new double[pMax - 1];
            for (int p = 0; p < pMax; ++p) {
                double[] ck = normAm(A, p+2);
                double c = Math.pow(ck[0], 1.0/(p+2));
                mv = mv + ck[1];
                eta[p] = c;
            }
            for (int p = 0; p < pMax - 1; ++p) {
                alpha[p] = Math.max(eta[p], eta[p+1]);
            }
        }
        DoubleMatrix M = DoubleMatrix.zeros(mMax, pMax-1);
        for (int p = 2; p <= pMax; ++p) {
            for (int m = p * (p-1) - 1; m <= mMax; ++m)
                M.put(m-1, p-2, alpha[p-2] / ThetaTaylor.THETA[m-1]);
        }
        return M;
    }
    
    private static final double[] normAm(final DoubleMatrix A, final int m) {
        int t = 1;
        final int n = A.getColumns();
        double c, mv;
        if (A.equals(A.copy().assign(Functions.abs))) {
            DoubleMatrix e = DoubleFactory2D.dense.make(n, 1, 1.0);
            for (int j = 0; j < m; ++j)
                e = Algebra.DEFAULT.mult(Algebra.DEFAULT.transpose(A), e);
            c = Algebra.DEFAULT.normInfinity(e);
            mv = m;
        } else {
            
            DoubleMatrixFunction afunPower =
                    new DoubleMatrixFunction() {
                        public DoubleMatrix apply(DoubleMatrix X,
                                boolean transpose) {
                            
                            if (!transpose) {
                                for (int i = 0; i < m; ++i) {
                                    X = Algebra.DEFAULT.mult(A, X);
                                }
                            } else {
                                DoubleMatrix AT =
                                        Algebra.DEFAULT.transpose(A);
                                for (int i = 0; i < m; ++i) {
                                    X = Algebra.DEFAULT.mult(AT, X);
                                } 
                            }
                            return X;
                        }
                        public int getDimensions() { return n; }
                        public boolean isReal() { return true; }
            };
            
            double[] c_it = HT00.normest1(afunPower, t);
            c = c_it[0];
            mv = c_it[1] * t * m;
        }
        return new double[]{c, mv};
    }
    
}
