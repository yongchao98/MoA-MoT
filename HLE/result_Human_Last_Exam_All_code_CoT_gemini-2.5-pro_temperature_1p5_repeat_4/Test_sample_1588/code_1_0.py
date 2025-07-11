import math

def get_prime_factorization(num):
    """Returns the prime factorization of num as a dictionary."""
    factors = {}
    d = 2
    temp_num = num
    while d * d <= temp_num:
        while temp_num % d == 0:
            factors[d] = factors.get(d, 0) + 1
            temp_num //= d
        d += 1
    if temp_num > 1:
        factors[temp_num] = factors.get(temp_num, 0) + 1
    return factors

def moebius(n):
    """Calculates the Moebius function mu(n)."""
    if n == 1:
        return 1
    factors = get_prime_factorization(n)
    for p in factors:
        if factors[p] > 1:
            return 0
    return (-1)**len(factors)

def get_divisors(n):
    """Returns a sorted list of all positive divisors of n."""
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n//i)
    return sorted(list(divs))

def solve_bch_coeffs():
    """
    Calculates the number of nonzero coefficients of a given order in the
    BCH expansion using Witt's formula.
    """
    n = 10
    k = 2

    print(f"Calculating the number of nonzero BCH coefficients for order n={n} and k={k} generators.")
    print(f"Using Witt's Formula: L_k(n) = (1/n) * Σ[d|n] (μ(d) * k^(n/d))")
    print("-" * 30)

    divisors = get_divisors(n)
    total_sum = 0
    
    # We will build the string for the final equation step by step
    equation_terms = []

    for d in divisors:
        mu_d = moebius(d)
        power_val = k**(n // d)
        term = mu_d * power_val
        total_sum += term
        
        sign = "+" if term >= 0 else "-"
        # For the first term, we omit the leading '+'
        if not equation_terms:
            if sign == "+":
                sign = ""
        else:
            sign = f" {sign} "
            
        # Add the detailed term to our equation string
        equation_terms.append(f"{sign}({mu_d}) * {k}^{n//d}")
        
    result = total_sum // n
    
    # Now build the numerical version of the equation
    numerical_terms = []
    current_sum = 0
    for d in divisors:
        mu_d = moebius(d)
        power_val = k**(n//d)
        term = mu_d * power_val
        current_sum += term
        
        sign = "+" if term >= 0 else "-"
        # For the first term, we omit the leading '+'
        if not numerical_terms:
             # handle case where first term is negative
            sign = "" if sign == "+" else "-"
        else:
            sign = f" {sign} "
        numerical_terms.append(f"{sign}{abs(term)}")
        
    print(f"The number of coefficients is (1/{n}) * ( " + " ".join(numerical_terms) + f" )")
    print(f"= (1/{n}) * ( {total_sum} )")
    print(f"= {result}")
    
    return result

if __name__ == "__main__":
    solve_bch_coeffs()
    print("\nThus, the number of nonzero coefficients of order 10 is 99.")
    print("<<<99>>>")
