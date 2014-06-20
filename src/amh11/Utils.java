/**
 * HT00.java
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

import java.util.Random;

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.MatrixEntry;
import no.uib.cipr.matrix.Vector;
import no.uib.cipr.matrix.VectorEntry;

/**
 * @author Arman D. Bilge <armanbilge@gmail.com>
 *
 */
public final class Utils {

    private Utils() {}
    
    public static interface MatrixFunction {
        public Matrix apply(Matrix X, boolean transpose);
        public int getDimensions();
        public boolean isReal();
    }
    
    public static interface AMH11Random {
        public boolean nextBoolean();
        public double nextDouble();
    }
    
    private static AMH11Random random = new AMH11Random() {
        private final Random random = new Random();
        public boolean nextBoolean() {
            return random.nextBoolean();
        }
        public double nextDouble() {
            return random.nextDouble();
        }
    };
    
    public static final void setRandom(AMH11Random random) {
        Utils.random = random;
    }
    
    public static final Vector fill(Vector v, double d) {
        for (int i = 0; i < v.size(); ++i) v.set(i, d);
        return v;
    }
    
    public static final double trace(Matrix M) {
        double trace = 0.0;
        for (int i = 0; i < M.numRows(); ++i)
            trace += M.get(i, i);
        return trace;
    }

    public static final Matrix ascendingDiagonal(int n) {
        Matrix M = new DenseMatrix(n, n);
        for (int i = 0; i < n; ++i) M.set(i, i, i+1);
        return M;
    }
    
    public static final Matrix ceil(Matrix M) {
        for (MatrixEntry me : M) me.set(Math.ceil(me.get()));
        return M;
    }
    
    public static final Matrix abs(Matrix M) {
        for (MatrixEntry me : M) me.set(Math.abs(me.get()));
        return M;
    }
    
    public static final Matrix Zero2Inf(Matrix M) {
        for (int i = 0; i < M.numRows(); ++i) {
            for (int j = 0; j < M.numColumns(); ++j)
                if (M.get(i, j) == 0) M.set(i, j, Double.POSITIVE_INFINITY);
        }
        return M;
    }
    
    public static final boolean isPositive(Matrix M) {
        for (MatrixEntry me : M) {
            if (me.get() < 0) return false;
        }
        return true;
    }
    
    public static final double entrySum(Vector v) {
        double sum = 0.0;
        for (VectorEntry ve : v) sum += ve.get();
        return sum;
    }
    
    public static final Vector rowAsVector(Matrix M, int row) {
        Vector v = new DenseVector(M.numColumns());
        for (int i = 0; i < M.numColumns(); ++i) v.set(i, M.get(row, i));
        return v;
    }

    public static final Vector columnAsVector(Matrix M, int column) {
        Vector v = new DenseVector(M.numRows());
        for (int i = 0; i < M.numColumns(); ++i) v.set(i, M.get(i, column));
        return v;
    }

    public static final Vector absColumnSums(Matrix M) {
        Vector sums = new DenseVector(M.numColumns());
        for (int i = 0; i < M.numColumns(); ++i) {
            double sum = Utils.columnAsVector(M, i).norm(Vector.Norm.One);
            sums.set(i, sum);
        }
        return sums;
    }
    
    public static final Matrix randomSigns(int r, int c) {
        Matrix signs = new DenseMatrix(r, c);
        for (int i = 0; i < r; ++i) {
            for (int j = 0; j < c; ++j) {
                double sign = random.nextBoolean() ? 1.0 : -1.0;
                signs.set(i, j, sign);
            }
        }
        return signs;
    }

    public static final Vector columnsMax(Matrix M) {
        Vector maxes = new DenseVector(M.numColumns());
        for (int i = 0; i < maxes.size(); ++i)
            maxes.set(i, Utils.columnAsVector(M, i).norm(Vector.Norm.Infinity));
        return maxes;
    }
    
    public static final int memberCount(int[] a, Vector b) {
        int count = 0;
        for (int i = 0; i < a.length; ++i) {
            double d = a[i];
            if (contains(b, d)) ++count;
        }
        return count;
    }
    
    public static final boolean contains(Vector v, double d) {
        for (int i = 0; i < v.size(); ++i) {
            if (v.get(i) == d) return true;
        }
        return false;
    }

    public static final Matrix createSubMatrix(Matrix M, int row, int column,
            int rows, int columns) {
        Matrix S = new DenseMatrix(rows, columns);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < columns; ++j)
                S.set(i, j, M.get(row+i, column+j));
        }
        return S;
    }
}
