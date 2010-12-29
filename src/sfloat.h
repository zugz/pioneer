#ifndef SFLOAT_H
#define SFLOAT_H

#include <SDL_stdinc.h>
#include <math.h>
#include <assert.h>

#ifndef CLAMP
# define CLAMP(a, min, max)      (((a) > (max)) ? (max) : (((a) < (min)) ? (min) : (a)))
#endif /* CLAMP */

#define SFLOAT_PVE_INF 0x7f800000
#define SFLOAT_NVE_INF 0xff800000

/**
 * A quick and dirty software floating point implementation. If you want something good look at SoftFloat
 * or some other decent implementation.
 */
class sfloat {
public:
	union {
		Uint32 val;
		float fval;
	};
	
	/* raw ieee float32 binary */
//	sfloat(Uint32 m): val(m) {}
//	sfloat(Sint32 m): val(m) {}
	sfloat(): val(0) {}
	/* XXX Be careful to avoid using float literals that can't be precisely
	   represented, like 0.1f for example */
//	sfloat(float v): fval(v) {}
	/** As fraction */
	sfloat(Uint64 numerator, Uint64 denominator) {
		Uint64 mantissa = (numerator<<32) / denominator;
		upck_float f;
		f.sign = 0;
		f.e = -32;
		while (mantissa >= (1ULL<<24)) {
			mantissa >>= 1;
			f.e++;
		}
		f.m = mantissa;
		f.Normalize();
		*this = Pack(f);
	}
	sfloat(unsigned int m, short e, bool sign) {
		upck_float f;
		f.m = m;
		f.e = e;
		f.sign = sign ? 1 : 0;
		f.Normalize();
		*this = Pack(f);
	}


	bool IsInfinite() const {
		return (val == 0xff800000) || (val == 0x7f800000);
	}
	bool IsNaN() const {
		return (((val >> 23) & 0xff) == 0xff) && (val & 0x7fffff);
	}

	friend bool operator==(sfloat a, sfloat b) {
		return a.val == b.val;
	}
	friend bool operator!=(sfloat a, sfloat b) { return a.val != b.val; }
	friend bool operator<(sfloat a, sfloat b) { return b>a; }
	friend bool operator>(sfloat a, sfloat b) {
		if (a == b) return false;
		return a >= b;
	}
	friend bool operator>=(sfloat _a, sfloat _b) {
		if (_a == _b) return true;
		if (_a.IsInfinite() || _b.IsInfinite()) return false;
		upck_float a = Unpack(_a);
		upck_float b = Unpack(_b);
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
		r.val = a.val ^ (1<<31);
		return r;
	}
	double ToDouble() const {
		return (double)ToFloat();
	}
	float ToFloat() const { return fval; }
	friend sfloat operator*(sfloat _a, sfloat _b) {
		if (_a.IsInfinite()) return _a;
		if (_b.IsInfinite()) return _b;
		if (_a.IsNaN()) return _a;
		if (_b.IsNaN()) return _b;
		upck_float a = Unpack(_a);
		upck_float b = Unpack(_b);
		upck_float r;
		r.e = a.e + b.e;
		Uint64 m = (Uint64)a.m * (Uint64)b.m;
		while (m >= (1LL<<24)) {
			m >>= 1;
			r.e++;
		}
		r.m = m;
		r.sign = (a.sign + b.sign) & 0x1;
		return Pack(r);
	}
	sfloat inverseOf() {
		return sfloat(1,0,false) / *this;
	}
	friend sfloat operator/(sfloat _a, sfloat _b) {
		upck_float a = Unpack(_a);
		upck_float b = Unpack(_b);
		upck_float r;
		if (b.m == 0) {
			sfloat r;
			r.val = SFLOAT_PVE_INF;
			return r;
		}

		r.e = a.e - 32 - b.e;
		Uint64 rm = (((Uint64)a.m)<<32) / (Uint64)b.m;
		while (rm >= (1ULL<<24)) {
			rm >>= 1;
			r.e++;
		}
		r.m = rm;
		r.sign = (a.sign + b.sign) & 0x1;
		return Pack(r);
	}
	void Print() const {
		//printf("DEBUG: %.15e (%s %u x2^ %hd)\n", ToDouble(), (sign ? "- " : "+ "), m, e);
		upck_float f = Unpack(*this);
		printf("sfloat(%uU,%hd,%s), // [0x%08x] %e\n", f.m, f.e, f.sign ? "true" : "false", val, ToDouble());
	}
	friend sfloat operator+(sfloat a, sfloat b) {
		return (a - (-b));
	}
	sfloat &operator-=(sfloat b) { *this = (*this) - b; return *this; }	
	sfloat &operator+=(sfloat b) { *this = (*this) + b; return *this; }	
	sfloat &operator*=(sfloat b) { *this = (*this) * b; return *this; }	
	sfloat &operator/=(sfloat b) { *this = (*this) / b; return *this; }	

	friend sfloat operator-(sfloat _a, sfloat _b) {
		if (_a.IsInfinite()) return _a;
		if (_b.IsInfinite()) return _b;
		if (_a.IsNaN()) return _a;
		if (_b.IsNaN()) return _b;
		upck_float a = Unpack(_a);
		upck_float b = Unpack(_b);
		if (!a.m) return -_b;
		if (!b.m) return _a;
		assert(Pack(a) == _a);
		assert(Pack(b) == _b);
		// lose some precision so they have matching exponents
		while (a.e < b.e) {
			a.m >>= 1;
			a.e++;
		}
		while (b.e < a.e) {
			b.m >>= 1;
			b.e++;
		}
		upck_float r;
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
		if (m >= (1LL<<24)) {
			m >>= 1;
			r.e++;
		}
		r.m = (unsigned int)m;
		r.Normalize();
		return Pack(r);
	}
	/** floor */
	int ToInt32() const {
		upck_float r = Unpack(*this);
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
		return r.sign ? -r.m : r.m;
	}
	Sint64 ToInt64() const { return ToInt32(); } // PISH
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
		sfloat r = *this;
		r.val &= 0x7fffffff;
		return r;
	}

private:
	static const sfloat oneThird;
	static const sfloat inverseFactorials[32];
	static const sfloat inverses[32];
	static const sfloat exp_pow2[8];
	static const sfloat exp_negpow2[7];
	
	struct upck_float {
		Uint32 m; // mantissa
		Sint16 e; // exponent
		Uint8 sign;
		void Normalize() {
			if (m) {
				while (m >= (1LL<<24)) {
					m >>= 1;
					e++;
				}
				while (!(m&(1<<23))) {
					m <<= 1;
					e--;
				}
			}
			if (!m) {
				e = 0;
			}
		}
	};

	static sfloat Pack(const upck_float &f) {
		sfloat r;
		if ((!f.m) || ((f.e+150) < 0)) {
			r.val = 0;
		} else if (f.e + 150 >= 255) {
			// infinity
			r.val = (((Uint32)f.sign)<<31) | (0xff << 23);
		} else {
			r.val = (((Uint32)f.sign)<<31) | (((f.e+150)&0xff)<<23) | (f.m & 0x7fffff);
		}
		return r;
	}
	static upck_float Unpack(sfloat v) {
		upck_float f;
		assert(!v.IsInfinite());
		assert(!v.IsNaN());
		if (!v.val) {
			f.m = f.e = f.sign = 0;
		} else {
			f.m = (v.val & 0x7fffff) | 0x800000;
			f.e = (signed)((v.val>>23)&0xff)-150;
			f.sign = (v.val>>31) & 1;
		}
		return f;
	}
};

#if 0
template<int highbitM, int restM, int E>
struct __sfloat {
		enum { mantissa = __sfloat< (restM>>30)&1, (restM<<1), E-1>::mantissa,
					exponent = __sfloat< (restM>>30)&1, (restM<<1), E-1>::exponent };
};
template<int restM, int E>
struct __sfloat<1, restM, E> {
		enum { mantissa = (1<<31) | restM, exponent = E };
};
#define static_sfloat(x) sfloat(__sfloat<(x>>31)&1,x,0>::mantissa, __sfloat<(x>>31)&1,x,0>::exponent, false)
#endif

#endif /* SFLOAT_H */
