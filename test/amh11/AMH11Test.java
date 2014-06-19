/**
 * AMH11Test.java
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

import static org.junit.Assert.assertTrue;

import java.util.Random;

import org.jblas.DoubleMatrix;
import org.jblas.MatrixFunctions;
import org.junit.Test;

import cern.colt.matrix.tdouble.DoubleFactory1D;
import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.jet.math.tdouble.DoubleFunctions;

public class AMH11Test {

    private static final Random random = new Random(123);
        
    @Test
    public void test() {
        
        int size = 8;
        for (int i = 0; i < 1024; ++i) {
            DoubleMatrix2D M = randomMatrix(size)
                    .assign(DoubleFunctions.mult(2.0))
                    .assign(DoubleFunctions.min(-1.0));
            DoubleMatrix1D v = randomVector(size);
            double t = Math.random();
            DoubleMatrix1D amh11 = AMH11.expmv(t, M, v);
            DoubleMatrix jblas = MatrixFunctions.expm(
                    new DoubleMatrix(M.toArray()).muli(t)).mmul(
                    new DoubleMatrix(v.toArray()));
            for (int j = 0; j < amh11.size(); ++j) {
                assertTrue(same(amh11.get(j), jblas.get(j)));
            }
        }
        
    }
    
    private static final DoubleMatrix2D randomMatrix(int size) {
        DoubleMatrix2D R = FlexibleDoubleFactory2D.large.make(size, size);
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j)
                R.set(i, j, random.nextDouble());
        }
        return R;
    }
    
    private static final DoubleMatrix1D randomVector(int size) {
        DoubleMatrix1D r = DoubleFactory1D.dense.make(size);
        for (int i = 0; i < size; ++i)
            r.set(i, random.nextDouble());
        return r;
    }
    
    @SuppressWarnings("unused")
    private static double EPSILON = 2.220446049250313E-16;
    private static double SQRT_EPSILON = 1.4901161193847656E-8;
    @SuppressWarnings("unused")
    private static double SQRT_SQRT_EPSILON = 1.220703125E-4;
    private static boolean same(double a, double b) {
        return Math.abs((a/b)-1.0) <= SQRT_EPSILON;
    }
    
}
