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
import org.apache.commons.math3.linear.RealMatrix;

public class AMH11 {

    private static final void normAm(RealMatrix A, double m) {
        int t = 1;
        int n = A.getColumnDimension();
        if (A.equals(abs(A))) {
            double[] ones = new double[n];
            Arrays.fill(ones, 1.0);
            RealMatrix e = new Array2DRowRealMatrix(ones);
            for (int j = 0; j < m; ++j)
                e = A.transpose().multiply(e);
            double c = e.getRowVector(0).getLInfNorm();
            double mv = m;
        } else {
            throw new RuntimeException("Not implemented!");
        }
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
