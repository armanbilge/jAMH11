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

import java.util.Arrays;
import java.util.Random;

import org.jblas.DoubleMatrix;
import org.jblas.MatrixFunctions;
import org.junit.Test;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;

public class AMH11Test {

    private final Random random = new Random();
    
    public AMH11Test() {
    }

    @Test
    public void test() {
        
        int size = 32;
        DoubleMatrix2D M = DoubleFactory2D.dense.random(size, size);
        DoubleMatrix1D v = DoubleFactory1D.dense.random(size);
        double t = random.nextDouble();
        
        DoubleMatrix1D amh11 = AMH11.expmv(t, M, v);
        DoubleMatrix jblas = MatrixFunctions.expm(new DoubleMatrix(M.toArray())
                .muli(t)).mmul(new DoubleMatrix(v.toArray()));
        
        System.out.println(Arrays.deepToString(M.toArray()));
        System.out.println(Arrays.toString(v.toArray()));
        System.out.println(t);
        
        System.out.println(Arrays.toString(amh11.toArray()));
        System.out.println(Arrays.toString(jblas.toArray()));
        for (int i = 0; i < amh11.size(); ++i) {
            assertTrue(same(amh11.getQuick(i), jblas.get(i)));
        }
    }
    
    private static double EPSILON = 2.220446049250313E-16;
    private static double SQRT_EPSILON = 1.4901161193847656E-8;
    private static double SQRT_SQRT_EPSILON = 1.220703125E-4;
    private static boolean same(double a, double b) {
        return Math.abs((a/b)-1.0) <= SQRT_SQRT_EPSILON;
    }
    
}
