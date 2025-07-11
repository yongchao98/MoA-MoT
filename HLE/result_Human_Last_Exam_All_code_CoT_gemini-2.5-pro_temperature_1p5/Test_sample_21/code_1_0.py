import sympy
from sympy import Poly, symbols, factor_list

def calculate_density():
    """
    This function numerically estimates the density of primes p for which
    the polynomial f(x) = x^7 - 14x^5 + 56x^3 - 56x + 22 is irreducible mod p.
    """
    # Set a limit for checking primes. A higher limit gives a better approximation.
    PRIME_LIMIT = 10000

    # Define the variable and the polynomial using sympy
    x = symbols('x')
    f = x**7 - 14*x**5 + 56*x**3 - 56*x + 22

    # Initialize counters
    # The Chebotarev Density Theorem applies to primes that do not ramify.
    # The ramified primes are the prime divisors of the discriminant (2, 7, 29).
    # We will exclude them from our count.
    prime_count = 0
    irreducible_count = 0

    print(f"Checking for primes up to {PRIME_LIMIT}...")
    
    # Iterate through primes up to the limit
    for p in sympy.primerange(2, PRIME_LIMIT):
        # Exclude ramified primes
        if p in [2, 7, 29]:
            continue
        
        prime_count += 1
        
        # Use sympy's factor_list to factor the polynomial modulo p.
        # It returns a tuple: (leading_coefficient, [(factor1, exp1), (factor2, exp2), ...])
        _coeff, factors = factor_list(f, modulus=p)

        # The polynomial is irreducible if there is only one factor with an exponent of 1.
        if len(factors) == 1 and factors[0][1] == 1:
            irreducible_count += 1

    # Calculate the numerical density
    numerical_density = irreducible_count / prime_count if prime_count > 0 else 0
    
    # The theoretical density from our analysis is 6/7
    theoretical_density = 6/7

    print(f"Total non-ramified primes checked: {prime_count}")
    print(f"Number of primes where f(x) is irreducible: {irreducible_count}")
    print(f"Numerical Density: {irreducible_count} / {prime_count} = {numerical_density:.6f}")
    print(f"Theoretical Density: 6 / 7 = {theoretical_density:.6f}")

if __name__ == '__main__':
    calculate_density()
