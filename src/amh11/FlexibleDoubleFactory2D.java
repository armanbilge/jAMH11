/**
 * FlexibleDoubleFactory2D.java
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

import cern.colt.matrix.tdouble.DoubleFactory2D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.impl.DenseLargeDoubleMatrix2D;

/**
 * @author Arman D. Bilge <armanbilge@gmail.com>
 *
 */
public class FlexibleDoubleFactory2D extends DoubleFactory2D {

    private static final long serialVersionUID = -662407369244803206L;

    public static final DoubleFactory2D large = new FlexibleDoubleFactory2D();
    
    protected FlexibleDoubleFactory2D() {}
    
    @Override
    public DoubleMatrix2D make(int rows, int columns) {
        try {
            return super.make(rows, columns);
        } catch (IllegalArgumentException exc) {
            if (!"matrix too large".equals(exc.getMessage()))
                throw exc;
            return new DenseLargeDoubleMatrix2D(rows, columns);
        }
    }
    
}
