#ifndef SFLOAT_H
#define SFLOAT_H

#include <SDL_stdinc.h>
#include <math.h>

class sfloat {
public:
	unsigned int m; // mantissa
	short e; // exponent
	unsigned char sign;
	unsigned char type;

	enum sfloat_type {
		TYPE_NUMBER = 0,
		TYPE_INF = 1,
		TYPE_NAN
	};

	/** As fraction */
	sfloat(Uint64 numerator, Uint64 denominator) {
		Uint64 mantissa = (numerator<<32) / denominator;
		sign = 0;
		e = -32;
		while (mantissa >= (1ULL<<32)) {
			mantissa >>= 1;
			e++;
		}
		m = mantissa;
		type = TYPE_NUMBER;
		Normalize();
	}

	/** Native */
	sfloat(unsigned int m, short e, bool sign, char type = TYPE_NUMBER) {
		this->m = m;
		this->e = e;
		this->sign = sign ? 1 : 0;
		this->type = type;
		Normalize();
	}
	sfloat(int m) {
		this->m = abs(m);
		this->e = 0;
		this->sign = (m>>31)&1;
		this->type = TYPE_NUMBER;
		Normalize();
	}
	sfloat() {
		m = e = sign = type = 0;
		Normalize();
	}
	friend bool operator==(sfloat a, sfloat b) {
		a.Normalize();
		b.Normalize();
		return ((a.type == b.type) && (a.m==b.m) && (a.e==b.e) && (a.sign==b.sign));
	}
	friend bool operator!=(sfloat a, sfloat b) { return !(a == b); }
	friend bool operator<(sfloat a, sfloat b) { return b>a; }
	friend bool operator>(sfloat a, sfloat b) {
		if (a == b) return false;
		return a >= b;
	}
	friend bool operator>=(sfloat a, sfloat b) {
		a.Normalize();
		b.Normalize();
		if (a == b) return true;
		if (a.sign) {
			if (!b.sign) return false;
			if (a.m == 0) return true;
			if (b.m == 0) return false;
			if (a.e > b.e) return false;
			if (a.e == b.e) {
				if (a.m > b.m) return false;
			}
			return true;
		} else {
			if (b.sign) return true;
			if (a.m == 0) return false;
			if (b.m == 0) return true;
			if (a.e > b.e) return true;
			if (a.e == b.e) {
				if (a.m > b.m) return true;
			}
			return false;
		}
	}
	friend bool operator<=(sfloat a, sfloat b) { return b >= a; }

	friend sfloat operator-(sfloat a) {
		sfloat r;
		r.m = a.m;
		r.e = a.e;
		r.sign = !a.sign;
		return r;
	}
	double ToDouble() const {
		if (type == TYPE_NUMBER) {
			if (sign) return -(double)m * pow(2.0, (double)e);
			else return (double)m * pow(2.0, (double)e);
		} else if (type == TYPE_INF) {
			return 1.0/0.0;
		} else {
			return 0.0/0.0; // NAN
		}
	}
	float ToFloat() const { return (float)ToDouble(); }
	int ToInt32() const { return floor().m; }
	Sint64 ToInt64() const { return floor().m; } // PISH
	friend sfloat operator*(sfloat a, sfloat b) {
		if ((a.type == TYPE_NAN) || (b.type == TYPE_NAN)) {
			return sfloat(0,0,false,TYPE_NAN);
		}
		sfloat r;
		r.e = a.e + b.e;
		Uint64 m = (Uint64)a.m * (Uint64)b.m;
		while (m >= (1LL<<32)) {
			m >>= 1;
			r.e++;
		}
		r.m = m;
		r.sign = (a.sign + b.sign) & 0x1;
		return r;
	}
	sfloat inverseOf() {
		return sfloat(1,0,false) / *this;
	}
	friend sfloat operator/(sfloat a, sfloat b) {
		if ((a.type == TYPE_NAN) || (b.type == TYPE_NAN)) {
			return sfloat(0,0,false,TYPE_NAN);
		}
		sfloat r;
		a.Normalize();
		b.Normalize();
		if (b.m == 0) {
			return sfloat(0,0,false,TYPE_INF);
		}

		r.e = a.e - 32 - b.e;
		Uint64 rm = (((Uint64)a.m)<<32) / (Uint64)b.m;
		while (rm >= (1ULL<<32)) {
			rm >>= 1;
			r.e++;
		}
		r.m = rm;
		r.sign = (a.sign + b.sign) & 0x1;
		return r;
	}
	void Normalize() {
		if (m) {
			while (!(m&(1<<31))) {
				m <<= 1;
				e--;
			}
		}
		if (!m) {
			e = -31;
		}
	}
	void Print() const {
		//printf("DEBUG: %.15e (%s %u x2^ %hd)\n", ToDouble(), (sign ? "- " : "+ "), m, e);
		printf("sfloat(%uU,%hd,%s), // %e\n", m, e, sign ? "true" : "false", ToDouble());
	}
	friend sfloat operator+(sfloat a, sfloat b) {
		b.sign = !b.sign;
		return (a - b);
	}
	friend sfloat operator+(sfloat a, int b) { return a + sfloat(b); }
	friend sfloat operator+(int a, sfloat b) { return b+a; }
	friend sfloat operator-(sfloat a, int b) { return a - sfloat(b); }
	friend sfloat operator-(int a, sfloat b) { return sfloat(a) - b; }
	friend sfloat operator*(sfloat a, int b) { return a * sfloat(b); }
	friend sfloat operator*(int a, sfloat b) { return b*a; }
	sfloat &operator+=(int b) { *this = (*this) + b; return *this; }
	sfloat &operator-=(int b) { *this = (*this) - b; return *this; }
	sfloat &operator*=(int b) { *this = (*this) * b; return *this; }
	sfloat &operator-=(sfloat b) { *this = (*this) - b; return *this; }	
	sfloat &operator+=(sfloat b) { *this = (*this) + b; return *this; }	
	sfloat &operator*=(sfloat b) { *this = (*this) * b; return *this; }	
	sfloat &operator/=(sfloat b) { *this = (*this) / b; return *this; }	

	friend sfloat operator-(sfloat a, sfloat b) {
		if ((a.type == TYPE_NAN) || (b.type == TYPE_NAN)) {
			return sfloat(0,0,false,TYPE_NAN);
		}
		if (!a.m) {
			b.sign = !b.sign;
			return b;
		}
		if (!b.m) return a;
		a.Normalize();
		b.Normalize();
		// lose some precision so they have matching exponents
		while (a.e < b.e) {
			a.m >>= 1;
			a.e++;
		}
		while (b.e < a.e) {
			b.m >>= 1;
			b.e++;
		}
		sfloat r;
		r.e = a.e;
		Uint64 m;
		if (a.m >= b.m) {
			if (a.sign != b.sign) {
				m = (Uint64)a.m + (Uint64)b.m;
			} else {
				m = (Uint64)a.m - (Uint64)b.m;
			}
			r.sign = a.sign;
		} else {
			if (a.sign != b.sign) {
				m = (Uint64)b.m + (Uint64)a.m;
			} else {
				m = (Uint64)b.m - (Uint64)a.m;
			}
			r.sign = !b.sign;
		}
		if (m >= (1LL<<32)) {
			m >>= 1;
			r.e++;
		}
		r.m = (unsigned int)m;
		r.Normalize();
		return r;
	}
	sfloat floor() const {
		// wipe off the fraction bits
		sfloat r = *this;
		if (r.e > 0) {
			while (r.e > 0) {
				r.m <<= 1;
				r.e--;
			}
		} else if (r.e < 0) {
			while (r.e < 0) {
				r.m >>= 1;
				r.e++;
			}
		}
		return r;
	}
	static sfloat Log(sfloat a);
	static sfloat Exp(sfloat val);
	static sfloat Pow(sfloat b, sfloat e) {
		return Exp(Log(b)*e);
	}
	static sfloat Sqrt(sfloat a) {
		return Pow(a, sfloat(1,-1, false));
	}
	static sfloat CubeRoot(sfloat a) {
		return Pow(a, oneThird);
	}
	sfloat Abs() const {
		return sfloat(m, e, false, type);
	}

private:
	static const sfloat oneThird;
	static const sfloat inverseFactorials[32];
	static const sfloat inverses[32];
	static const sfloat exp_pow2[11];
	static const sfloat exp_negpow2[10];
};

#endif /* SFLOAT_H */
