from fractions import Fraction
import math

class PolyH:
    """A class for polynomial arithmetic in h, modulo h^5."""
    def __init__(self, coeffs):
        # coeffs are [c0, c1, c2, c3, c4] for c0 + c1*h + c2*h^2 + ...
        self.coeffs = [Fraction(c) for c in coeffs]
        if len(self.coeffs) > 5:
            self.coeffs = self.coeffs[:5]
        while len(self.coeffs) < 5:
            self.coeffs.append(Fraction(0))

    def __add__(self, other):
        return PolyH([self.coeffs[i] + other.coeffs[i] for i in range(5)])

    def __sub__(self, other):
        return PolyH([self.coeffs[i] - other.coeffs[i] for i in range(5)])

    def __mul__(self, other):
        new_coeffs = [Fraction(0)] * 5
        for i in range(5):
            for j in range(5):
                if i + j < 5:
                    new_coeffs[i + j] += self.coeffs[i] * other.coeffs[j]
        return PolyH(new_coeffs)

    def __rtruediv__(self, other):
        return PolyH([other]) / self

    def __truediv__(self, other):
        return self * other.inv()
        
    def __pow__(self, n):
        if n == 0:
            return PolyH([1])
        if n < 0:
            return self.inv() ** (-n)
        res = PolyH([1])
        base = self
        while n > 0:
            if n % 2 == 1:
                res = res * base
            base = base * base
            n //= 2
        return res

    def inv(self):
        # Assumes self.coeffs[0] == 1
        if self.coeffs[0] != 1:
            raise ValueError("Can only invert polynomials with constant term 1")
        
        # (1+x)^-1 = 1 - x + x^2 - x^3 + x^4 ...
        x = PolyH([0] + self.coeffs[1:])
        x_pow = PolyH([1])
        inv_poly = PolyH([1])
        for i in range(1, 5):
            x_pow = x_pow * x
            if i % 2 == 1:
                inv_poly = inv_poly - x_pow
            else:
                inv_poly = inv_poly + x_pow
        return inv_poly
        
    def __repr__(self):
        return f"PolyH({self.coeffs})"

    def to_string(self):
        parts = []
        if self.coeffs[0] != 0 or len(parts) == 0:
             parts.append(f"{self.coeffs[0]}")
        if self.coeffs[1] != 0:
             parts.append(f"{self.coeffs[1]}*h")
        for i in range(2, 5):
             if self.coeffs[i] != 0:
                parts.append(f"{self.coeffs[i]}*h^{i}")
        return " + ".join(parts).replace("+ -", "- ")

def nCr_f(n, r):
    return math.factorial(n) // math.factorial(r) // math.factorial(n-r)

def c_to_ch(c_poly, rank):
    """Converts a total Chern class polynomial to a Chern character polynomial."""
    c = c_poly.coeffs
    c1h = PolyH([0, c[1]])
    c2h2 = PolyH([0, 0, c[2]])
    c3h3 = PolyH([0, 0, 0, c[3]])
    c4h4 = PolyH([0, 0, 0, 0, c[4]])
    
    ch0 = PolyH([rank])
    ch1 = c1h
    ch2 = Fraction(1, 2) * (c1h * c1h - 2 * c2h2)
    ch3 = Fraction(1, 6) * (c1h*c1h*c1h - 3 * c1h * c2h2 + 3 * c3h3)
    ch4 = Fraction(1, 24) * (c1h**4 - 4*(c1h**2)*c2h2 + 2*(c2h2**2) + 4*c1h*c3h3 - 4*c4h4)
    
    return ch0 + ch1 + ch2 + ch3 + ch4

def ch_to_c(ch_poly):
    """Converts a Chern character polynomial back to a total Chern class."""
    ch = ch_poly.coeffs
    p1h = PolyH([0, ch[1]])
    p2h2 = PolyH([0, 0, 2 * ch[2]])
    p3h3 = PolyH([0, 0, 0, 6 * ch[3]])
    p4h4 = PolyH([0, 0, 0, 0, 24 * ch[4]])
    
    c1h = p1h
    c2h2 = Fraction(1, 2) * (c1h * p1h - p2h2)
    c3h3 = Fraction(1, 3) * (c1h * p2h2 - c2h2 * p1h + p3h3)
    c4h4 = Fraction(1, 4) * (c1h * p3h3 - c2h2 * p2h2 + c3h3 * p1h - p4h4)
    
    return PolyH([1]) + c1h + c2h2 + c3h3 + c4h4

def solve():
    """Main function to perform the calculation."""
    h = PolyH([0, 1])
    one = PolyH([1])

    # 1. Chern classes of T, Lambda^2 T, Lambda^3 T
    c_T = (one + h)**5
    c_L2T = (one + 2*h)**nCr_f(5,2) / c_T
    c_L3T = (one + 3*h)**nCr_f(5,3) / c_L2T

    # 2. Convert to Chern characters
    ch_T = c_to_ch(c_T, rank=4)
    ch_L2T = c_to_ch(c_L2T, rank=nCr_f(4,2))
    ch_L3T = c_to_ch(c_L3T, rank=nCr_f(4,3))

    # 3. Compute ch(S^(2,1)T)
    ch_S21T = ch_T * ch_L2T - ch_L3T
    
    # 4. Convert back to Chern class
    c_S21T = ch_to_c(ch_S21T)
    
    coeffs = [c.numerator if c.denominator == 1 else f"{c.numerator}/{c.denominator}" for c in c_S21T.coeffs]
    
    result = f"c(S^(2,1)T) = {coeffs[0]}"
    if coeffs[1] != 0:
        result += f" + {coeffs[1]}*h"
    if coeffs[2] != 0:
        result += f" + {coeffs[2]}*h^2"
    if coeffs[3] != 0:
        result += f" + {coeffs[3]}*h^3"
    if coeffs[4] != 0:
        result += f" + {coeffs[4]}*h^4"
    
    print(result.replace('+ -', '- '))

solve()
