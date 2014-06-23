package amh11;

import cern.colt.matrix.tdouble.DoubleFactory2D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.impl.DenseLargeDoubleMatrix2D;

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
