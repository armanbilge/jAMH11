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

import cern.colt.function.tdouble.DoubleFunction;
import cern.colt.matrix.tdouble.DoubleFactory1D;
import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.algo.DenseDoubleAlgebra;
import cern.jet.math.tdouble.DoubleFunctions;

/**
 * @author Arman D. Bilge <armanbilge@gmail.com>
 *
 */
public final class AMH11 {
    
    private static final double TOL = Math.pow(2, -53);
    
    private AMH11() {}
    
    public static final double[] expmv(double t, double[][] A, double[] b) {
        return expmv(t, FlexibleDoubleFactory2D.large.make(A),
                DoubleFactory1D.dense.make(b), null, true, false, true, false)
                .toArray();
    }
    
    public static final DoubleMatrix1D expmv(double t, DoubleMatrix2D A,
            DoubleMatrix1D b) {
        return expmv(t, A, b, null, true, false, true, true);
    }
    
    public static final DoubleMatrix1D expmv(double t, DoubleMatrix2D A,
            DoubleMatrix1D b, DoubleMatrix2D M, boolean shift, boolean bal,
            boolean fullTerm, boolean copy) {
           
        if (copy) {
            A = A.copy();
            b = b.copy();
        }
        
        if (bal) {
            throw new RuntimeException("Not implemented!");
        }
        
        int n = A.rows();
        double mu = 0.0;
        if (shift) {
            mu = DenseDoubleAlgebra.DEFAULT.trace(A) / n;
            A.assign(FlexibleDoubleFactory2D.large.identity(n)
                    .assign(DoubleFunctions.mult(mu)), DoubleFunctions.minus);
        }
        double tt;
        if (M == null) {
            tt = 1.0;
            M = selectTaylorDegree(
                    A.copy().assign(DoubleFunctions.mult(t)),
                    b, 55, 8, shift, bal, false);
        } else {
            tt = t;
        }
        
        double s = 1.0;
        int m, mMax;
        if (t == 0) {
            m = 0;
        } else {
            mMax = M.rows();
            DoubleMatrix2D U = FlexibleDoubleFactory2D.large
                    .diagonal(DoubleFactory1D.dense.ascending(mMax));
            DoubleMatrix2D C = DenseDoubleAlgebra.DEFAULT.mult(DenseDoubleAlgebra.DEFAULT
                    .transpose(M.assign(DoubleFunctions.mult(Math.abs(tt)))
                            .assign(DoubleFunctions.ceil)), U);
            C.assign(ZERO2Inf);
            double cost;
            int[] min = new int[2];
            cost = getMin(C, min);
            m = min[C.columns() > 1 ? 1 : 0] + 1;
            if (cost == Double.POSITIVE_INFINITY)
                cost = 0;
            s = Math.max(cost/m, 1);
        }
        
        double eta = 1;
        if (shift) eta = Math.exp(t * mu / s);
        DoubleMatrix1D f = b.copy();
        for (int i = 0; i < s; ++i) {
            double c1 = DenseDoubleAlgebra.DEFAULT.normInfinity(b);
            for (int k = 1; k <= m; ++k) {
                b = DenseDoubleAlgebra.DEFAULT.mult(A, b).assign(DoubleFunctions.mult(t/(s*k)));
                f.assign(b, DoubleFunctions.plus);
                double c2 = DenseDoubleAlgebra.DEFAULT.normInfinity(b);
                if (!fullTerm) {
                    if (c1 + c2 <= TOL * DenseDoubleAlgebra.DEFAULT.normInfinity(f))
                        break;
                    c1 = c2;
                }
            }
            b = f.assign(DoubleFunctions.mult(eta));
        }
        return f;
    }
    
    private static final DoubleFunction ZERO2Inf = new DoubleFunction() {
        public double apply(double d) {
            if (d == 0) return Double.POSITIVE_INFINITY;
            return d;
        }
    };
    
    private static final double getMin(DoubleMatrix2D M, int[] index) {
        double min = Double.POSITIVE_INFINITY;
        for (int i = 0; i < M.rows(); ++i) {
            for (int j = 0; j < M.columns(); ++j) {
                double d = M.getQuick(i, j);
                if (d < min) {
                    min = d;
                    index[0] = i;
                    index[1] = j;
                }
            }
        }
        return min;
    }
    
    public static final DoubleMatrix2D selectTaylorDegree(DoubleMatrix2D A,
            DoubleMatrix1D b, int mMax, int pMax, boolean shift,
            boolean bal, boolean forceEstm) {
        int n = A.rows();
        if (bal) {
            throw new RuntimeException("Not implemented!");
        }
        double mu;
        if (shift) {
            mu = DenseDoubleAlgebra.DEFAULT.trace(A) / n;
            A.assign(FlexibleDoubleFactory2D.large.identity(n)
                    .assign(DoubleFunctions.mult(-mu)), DoubleFunctions.plus);
        }
        double mv = 0.0;
        double normA = 0.0;
        if (!forceEstm)
            normA = DenseDoubleAlgebra.DEFAULT.norm1(A);
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
        DoubleMatrix2D M = FlexibleDoubleFactory2D.large.make(mMax, pMax-1, 0.0);
        for (int p = 2; p <= pMax; ++p) {
            for (int m = p * (p-1) - 1; m <= mMax; ++m)
                M.setQuick(m-1, p-2, alpha[p-2] / ThetaTaylor.THETA[m-1]);
        }
        return M;
    }
    
    private static final double[] normAm(final DoubleMatrix2D A, final int m) {
//        int t = 1;
        final int n = A.columns();
        double c, mv;
//        if (A.equals(A.copy().assign(DoubleFunctions.abs))) {
            DoubleMatrix2D e = FlexibleDoubleFactory2D.large.make(n, 1, 1.0);
            for (int j = 0; j < m; ++j)
                e = DenseDoubleAlgebra.DEFAULT.mult(DenseDoubleAlgebra.DEFAULT.transpose(A), e);
            c = DenseDoubleAlgebra.DEFAULT.normInfinity(e);
            mv = m;
//        } else {
//            
//            DoubleMatrix2DFunction afunPower =
//                    new DoubleMatrix2DFunction() {
//                        public DoubleMatrix2D apply(DoubleMatrix2D X,
//                                boolean transpose) {
//                            
//                            if (!transpose) {
//                                for (int i = 0; i < m; ++i) {
//                                    X = DenseDoubleAlgebra.DEFAULT.mult(A, X);
//                                }
//                            } else {
//                                DoubleMatrix2D AT =
//                                        DenseDoubleAlgebra.DEFAULT.transpose(A);
//                                for (int i = 0; i < m; ++i) {
//                                    X = DenseDoubleAlgebra.DEFAULT.mult(AT, X);
//                                } 
//                            }
//                            return X;
//                        }
//                        public int getDimensions() { return n; }
//                        public boolean isReal() { return true; }
//            };
//            
//            double[] c_it = HT00.normest1(afunPower, t);
//            c = c_it[0];
//            mv = c_it[1] * t * m;
//        }
        return new double[]{c, mv};
    }
    
}
