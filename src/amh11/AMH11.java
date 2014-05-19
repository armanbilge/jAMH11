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

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;

public class AMH11 {

    
    
    private static final Object[] selectTaylorDegree(RealMatrix A, RealMatrix b, 
            int mMax, int pMax, double prec, boolean shift, boolean bal,
            boolean forceEstm) {
        int n = A.getRowDimension();
        if (bal) {
            throw new RuntimeException("Not implemented!");
        }
        double mu;
        if (shift) {
            mu = A.getTrace() / n;
            A = A.add(MatrixUtils.createRealIdentityMatrix(n)
                    .scalarMultiply(-mu));
        }
        double mv = 0.0;
        double normA = 0.0;
        if (!forceEstm)
            normA = A.getNorm();
        int unA;
        double[] alpha;
        if (!forceEstm && normA < 4 * ThetaTaylor.THETA[mMax] * pMax
                * (pMax + 3) / (mMax * b.getColumnDimension())) {
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
        RealMatrix M = new Array2DRowRealMatrix(mMax, pMax-1);
        for (int p = 1; p < pMax; ++p) {
            for (int m = p * (p-1) - 1; m < mMax; ++m)
                M.setEntry(m, p-1, alpha[p-1] / ThetaTaylor.THETA[m]);
        }
        return new Object[]{M, mv, alpha, unA};
    }
    
    private static final double[] normAm(RealMatrix A, int m) {
        int t = 1;
        int n = A.getRowDimension();
        double c, mv;
        if (A.equals(abs(A))) {
            double[] ones = new double[n];
            Arrays.fill(ones, 1.0);
            RealMatrix e = new Array2DRowRealMatrix(ones);
            for (int j = 0; j < m; ++j)
                e = A.transpose().multiply(e);
            c = e.getRowVector(0).getLInfNorm();
            mv = m;
        } else {
            throw new RuntimeException("Not implemented!");
        }
        return new double[]{c, mv};
    }
        
    private static final RealMatrix abs(RealMatrix A) {
        RealMatrix B = A.copy();
        for (int i = 0; i < B.getRowDimension(); ++i) {
            for (int j = 0; j < B.getColumnDimension(); ++j)
                B.setEntry(i, j, Math.abs(B.getEntry(i, j)));
        }
        return B;
    }
    
}
