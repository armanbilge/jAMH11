/**
 * ThetaTaylor.java
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

import org.jblas.DoubleMatrix;

/**
 * @author Arman D. Bilge <armanbilge@gmail.com>
 *
 */
public final class Utils {

    private Utils() {}

    public static final double trace(DoubleMatrix A) {
        double trace = 0.0;
        for (int i = 0; i < A.getRows(); ++i) trace += A.get(i, i);
        return trace;
    }
    
    public static final DoubleMatrix ascendingMatrix(int n) {
        DoubleMatrix A = DoubleMatrix.zeros(n, n);
        for (int i = 0; i < n; ++i) A.put(i, i, i+1);
        return A;
    }
    
    public static final DoubleMatrix ceil(DoubleMatrix A) {
        for (int i = 0; i < A.getRows(); ++i) {
            for (int j = 0; j < A.getColumns(); ++j)
                A.put(i, j, Math.ceil(A.get(i, j)));
        }
        return A;
    }
    
    public static final DoubleMatrix abs(DoubleMatrix A) {
        for (int i = 0; i < A.getRows(); ++i) {
            for (int j = 0; j < A.getColumns(); ++j)
                A.put(i, j, Math.abs(A.get(i, j)));
        }
        return A;
    }
    
    public static final DoubleMatrix Zero2Inf(DoubleMatrix A) {
        for (int i = 0; i < A.getRows(); ++i) {
            for (int j = 0; j < A.getColumns(); ++j) {
                if (A.get(i, j) == 0.0) A.put(i, j, Double.POSITIVE_INFINITY);
            }
        }
        return A;        
    }
    
}
