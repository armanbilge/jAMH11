/**
 * AMH11Test.java
 *
 * AMH11: Java implementation of the matrix exponential method
 *     described by Al-Mohy and Higham (2011)
 *
 * Copyright (c) 2014 Arman Bilge <armanbilge@gmail.com>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

package amh11;

import static org.junit.Assert.assertTrue;

import java.util.Random;

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Matrices;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;

import org.jblas.DoubleMatrix;
import org.jblas.MatrixFunctions;
import org.junit.Test;

public class AMH11Test {

    public AMH11Test() {}

    private static final Random random = new Random(123);

    @Test
    public void test() {

        int size = 16;
        Matrix O = new DenseMatrix(size, size);
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) O.set(j, i, 1.0);
        }

        for (int i = 0; i < 256; ++i) {
            Matrix M = randomMatrix(size).scale(2).add(-1, O);
            Vector v = randomVector(size);
            double t = random.nextDouble();
            Vector amh11 = AMH11.expmv(t, M, v);
            DoubleMatrix jblas = MatrixFunctions.expm(
                    new DoubleMatrix(Matrices.getArray(M)).muli(t)).mmul(
                    new DoubleMatrix(Matrices.getArray(v)));
            for (int j = 0; j < amh11.size(); ++j) {
                assertTrue(same(amh11.get(j), jblas.get(j)));
            }
        }

    }

    private static final Matrix randomMatrix(int size) {
        Matrix R = new DenseMatrix(size, size);
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j)
                R.set(i, j, random.nextDouble());
        }
        return R;
    }

    private static final Vector randomVector(int size) {
        Vector r = new DenseVector(size);
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
