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

import no.uib.cipr.matrix.Matrices;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;

import org.jblas.DoubleMatrix;
import org.jblas.MatrixFunctions;
import org.junit.Test;

public class AMH11Test {

    private final Random random = new Random();
    
    public AMH11Test() {
    }

    @Test
    public void test() {
        
        int size = 16;
        Matrix M = Matrices.random(size, size);
        Vector v = Matrices.random(size);
        double t = random.nextDouble();
        
        Vector amh11 = AMH11.expmv(t, M, v);

        DoubleMatrix jblas = MatrixFunctions.expm(new DoubleMatrix(Matrices.getArray(M))
                .muli(t)).mmul(new DoubleMatrix(Matrices.getArray(v)));
        
        System.out.println(Arrays.deepToString(Matrices.getArray(M)));
        System.out.println(Arrays.toString(Matrices.getArray(v)));
        System.out.println(t);
        
        System.out.println(Arrays.toString(Matrices.getArray(amh11)));
        System.out.println(Arrays.toString(jblas.toArray()));
        for (int i = 0; i < amh11.size(); ++i) {
            assertTrue(same(amh11.get(i), jblas.get(i)));
        }
    }
    
    private static double EPSILON = 2.220446049250313E-16;
    private static double SQRT_EPSILON = 1.4901161193847656E-8;
    private static double SQRT_SQRT_EPSILON = 1.220703125E-4;
    private static boolean same(double a, double b) {
        return Math.abs((a/b)-1.0) <= SQRT_SQRT_EPSILON;
    }
    
}
