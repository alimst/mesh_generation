// User defined class named "rational" for rational numbers.
//#include <iostream>
#include <cmath>

class rational
{
    public:
	// construct a rational number from two int numbers, first as numerator and second as denominator.
            rational(int = 0 , int = 1); 

	// addition between two rational numbers.
            rational operator+(const rational& ) const;

	// addition beween an integer as first argument and a rational number as second agument.
	    friend rational operator+(const int& , const rational& );
 
	// subtraction between two rational numbers.
            rational operator-(const rational& ) const;

	// multiplication between two rational numbers.
            rational operator*(const rational& ) const;

	// multiplication between an integer as first argument and a rational number as second agument.
	    friend rational operator*(const int& , const rational& );

	// division between two rational numbers.	
            rational operator/(const rational& ) const;

	// less than or equal operation.
	    const bool operator<=(const rational& ) const;

	// assignment between two rational numbers.
	    rational& operator=(const rational& );
	// assignment of a rational number to an integer.
	    rational& operator=(const int& );

	// pre-increment operator
	    rational& operator++();

	// ceil function of a rational number.
	    static int ceil(const rational& );

	// floor function of a rational number.
	    static int floor(const rational& );

	// standard output operation
	    friend std::ostream& operator << (std::ostream&, const rational&);

	// standard input operation.
	    friend std::istream& operator >> (std::istream&, rational&);
 
    private:
            int num; // numerator
	        int den; // denominator
};

rational::rational(int numerator , int denominator )
{
	num = numerator;
	den = denominator;	
}

rational rational::operator+(const rational& y) const
{
	int resultNum = (y.num*den)+(y.den*num);
	int resultDen = y.den*den;
	return rational(resultNum , resultDen);
}

rational operator+(const int& n, const rational& y)
{
	rational x(n);
	return (x+y);
}

rational rational::operator-(const rational& y) const
{
	int resultNum = (y.den*num)-(y.num*den);
	int resultDen = y.den*den;
	return rational(resultNum , resultDen);
}
rational rational::operator*(const rational& y) const
{
	int resultNum = y.num*num;
	int resultDen = y.den*den;
	return rational(resultNum , resultDen);
}

rational operator*(const int& n, const rational& y)
{
	rational x(n);
	return (x*y);
}

rational rational::operator/(const rational& y) const
{
	int resultNum = num*y.den;
	int resultDen = den*y.num;
	return rational(resultNum , resultDen);
}

const bool rational::operator<=(const rational& y) const
{
	return ( (y.num*den) <= (y.den*num) );
}

rational& rational::operator=(const rational& y)
{
	if ( this != &y )
	{
		num = y.num;
		den = y.den;
	}
	return *this;
}

rational& rational::operator=(const int& y)
{
	num = y;
	den = 1;
	return *this;
}

rational& rational::operator++()
{
	num = num + den;
	return *this;
}

int rational::ceil(const rational& y) // implemented for only non-negative numbers
{
	if ( fmod(y.num , y.den) == 0 )
	{
		return ( y.num/y.den );
	}
	else
	{
		return ( (y.num+y.den-fmod(y.num , y.den))/y.den );
	}
}

int rational::floor(const rational& y) // implemented for only non-negative numbers
{
	return ( (y.num - fmod(y.num,y.den)) / y.den );
}

std::ostream& operator << (std::ostream& out, const rational& y)
{
	return out << y.num << "/" << y.den;
}

std::istream& operator >> (std::istream& in, rational& y)
{
	in >> y.num >> y.den;
	return in;
}
