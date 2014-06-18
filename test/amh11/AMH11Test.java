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

import java.lang.reflect.Field;
import java.util.Random;

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.Matrices;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;

import org.jblas.DoubleMatrix;
import org.jblas.MatrixFunctions;
import org.junit.BeforeClass;
import org.junit.Test;

public class AMH11Test {
    
    public AMH11Test() {}

    @BeforeClass
    public static void setup() throws Exception {
        Field field = Math.class.getDeclaredField("randomNumberGenerator");
        field.setAccessible(true);
        field.set(null, new Random(123));
    }
    
    @Test
    public void test() {
                
        int size = 16;
        Matrix O = new DenseMatrix(size, size);
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) O.set(j, i, 1.0);
        }
        
        for (int i = 0; i < 256; ++i) {
            Matrix M = Matrices.random(size, size).scale(2).add(-1, O);
            Vector v = Matrices.random(size);
            double t = Math.random();
            Vector amh11 = AMH11.expmv(t, M, v);
            DoubleMatrix jblas = MatrixFunctions.expm(
                    new DoubleMatrix(Matrices.getArray(M)).muli(t)).mmul(
                    new DoubleMatrix(Matrices.getArray(v)));
            for (int j = 0; j < amh11.size(); ++j) {
                assertTrue(same(amh11.get(j), jblas.get(j)));
            }
        }
        
    }
    
    private static double EPSILON = 2.220446049250313E-16;
    private static double SQRT_EPSILON = 1.4901161193847656E-8;
    private static double SQRT_SQRT_EPSILON = 1.220703125E-4;
    private static boolean same(double a, double b) {
        return Math.abs((a/b)-1.0) <= SQRT_EPSILON;
    }
    
}
