/**
 * AMH11.java
 * 
 * AMH11: Java implementation of Al-Mohy and Higham's (2011) matrix exponential
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

import cern.colt.function.DoubleFunction;
import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.linalg.Algebra;
import cern.colt.matrix.linalg.SmpBlas;
import cern.jet.math.Functions;

public class AMH11 {
    
    private static final double TOL = Math.pow(2, -53);
    
    public static final DoubleMatrix1D expmv(int t, DoubleMatrix2D A, DoubleMatrix1D b) {
        return expmv(t, A, b, null, true, false, true);
    }
    
    public static final DoubleMatrix1D expmv(int t, DoubleMatrix2D A,
            DoubleMatrix1D b, DoubleMatrix2D M, boolean shift, boolean bal,
            boolean fullTerm) {
        
        if (bal) {
            throw new RuntimeException("Not implemented!");
        }
        
        int n = A.columns();
        double mu = 0.0;
        if (shift) {
            mu = Algebra.DEFAULT.trace(A) / n;
            A.assign(DoubleFactory2D.rowCompressed.identity(n)
                    .assign(Functions.mult(-mu)), Functions.plus);
        }
        
        double tt, mv, mvd;
        if (M == null) {
            tt = 1.0;
            // TODO!
            Object[] M_mvd_alpha_unA = selectTaylorDegree(
                    A.assign(Functions.mult((double) t)),
                    b, Integer.MIN_VALUE, Integer.MIN_VALUE, shift, bal, false);
            M = (DoubleMatrix2D) M_mvd_alpha_unA[0];
            mv = (double) M_mvd_alpha_unA[1];
        } else {
            tt = t;
            mv = 0;
            mvd = 0;
        }
        
        double s = 1.0;
        int m, mMax, p;
        if (t == 0) {
            m = 0;
        } else {
            mMax = M.columns();
            p = M.rows();
            DoubleMatrix2D U = DoubleFactory2D.rowCompressed
                    .diagonal(DoubleFactory1D.dense.ascending(mMax));
            DoubleMatrix2D C = Algebra.DEFAULT.mult(Algebra.DEFAULT
                    .transpose(M.assign(Functions.mult(Math.abs(tt)))
                            .assign(Functions.abs)),  U);
            C.assign(ZERO2Inf);
            double cost;
            int[] min = new int[2];
            cost = getMin(C, min);
            m = min[C.columns() > 1 ? 1 : 0];
            if (cost == Double.POSITIVE_INFINITY)
                cost = 0;
            s = Math.max(cost/m, 1);
        }
        
        double eta = 1;
        if (shift) eta = Math.exp(t * mu /s);
        DoubleMatrix1D f = b.copy();
        DoubleMatrix1D intermediate = new DenseDoubleMatrix1D(f.size());
        for (int i = 0; i < s; ++i) {
            double c1 = Algebra.DEFAULT.normInfinity(b);
            for (int k = 1; k <= m; ++k) {
                SmpBlas.smpBlas.dgemv(false, t/(s*k), A, b, 0, intermediate);
                b = intermediate;
                mv += 1;
                f.assign(b, Functions.plus);
                double c2 = Algebra.DEFAULT.normInfinity(b);
                if (!fullTerm) {
                    if (c1 + c2 <= TOL * Algebra.DEFAULT.normInfinity(f))
                        break;
                    c1 = c2;
                }
            }
            b = f.assign(Functions.mult(eta));
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
        for (int i = 0; i < M.columns(); ++i) {
            for (int j = 0; j < M.rows(); ++j) {
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
    
    public static final Object[] selectTaylorDegree(DoubleMatrix2D A,
            DoubleMatrix1D b, int mMax, int pMax, boolean shift,
            boolean bal, boolean forceEstm) {
        int n = A.columns();
        if (bal) {
            throw new RuntimeException("Not implemented!");
        }
        double mu;
        if (shift) {
            mu = Algebra.DEFAULT.trace(A) / n;
            A.assign(DoubleFactory2D.rowCompressed.identity(n)
                    .assign(Functions.mult(-mu)), Functions.plus);
        }
        double mv = 0.0;
        double normA = 0.0;
        if (!forceEstm)
            normA = Algebra.DEFAULT.norm1(A);
        int unA;
        double[] alpha;
        if (!forceEstm && normA < 4 * ThetaTaylor.THETA[mMax] * pMax
                * (pMax + 3) / (mMax * b.size())) {
            unA = 1;
            double c = normA;
            alpha = new double[pMax - 1];
            Arrays.fill(alpha, 1.0);
        } else {
            unA = 0;
            double[] eta = new double[pMax - 1];
            alpha = new double[pMax - 1];
            for (int p = 0; p < pMax; ++p) {
                double[] ck = normAm(A, p+1);
                double c = Math.pow(ck[0], 1/(p+1));
                mv = mv + ck[1];
                eta[p] = c;
            }
            for (int p = 0; p < pMax - 1; ++p) {
                alpha[p] = Math.max(eta[p], eta[p+1]);
            }
        }
        DoubleMatrix2D M = DoubleFactory2D.dense.make(mMax, pMax-1, 0.0);
        for (int p = 1; p < pMax; ++p) {
            for (int m = p * (p-1) - 1; m < mMax; ++m)
                M.setQuick(m, p-1, alpha[p-1] / ThetaTaylor.THETA[m]);
        }
        return new Object[]{M, mv, alpha, unA};
    }
    
    private static final double[] normAm(DoubleMatrix2D A, int m) {
        int t = 1;
        int n = A.columns();
        double c, mv;
        if (A.equals(A.copy().assign(Functions.abs))) {
            DoubleMatrix2D e = DoubleFactory2D.dense.make(n, 1, 1.0);
            for (int j = 0; j < m; ++j)
                e = Algebra.DEFAULT.mult(Algebra.DEFAULT.transpose(A), e);
            c = Algebra.DEFAULT.normInfinity(e);
            mv = m;
        } else {
            throw new RuntimeException("Not implemented!");
        }
        return new double[]{c, mv};
    }
    
}
