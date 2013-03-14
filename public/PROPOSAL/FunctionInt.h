/**
 * function to be interpolated must be included in a class that implements this interface
 */
#ifndef FunctionInt_H
#define FunctionInt_H

class FunctionInt
{

public:

    virtual double functionInt(double x) const = 0;

//----------------------------------------------------------------------------//

    virtual ~FunctionInt(){}
    		

};

#endif //FunctionInt_H
