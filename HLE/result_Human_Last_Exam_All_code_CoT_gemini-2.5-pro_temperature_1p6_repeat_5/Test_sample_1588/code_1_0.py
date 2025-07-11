import math

def moebius(m):
    """
    Calculates the Möbius function mu(m).
    - mu(1) = 1
    - mu(m) = 0 if m has a squared prime factor.
    - mu(m) = (-1)^p if m is a product of p distinct primes.
    """
    if m == 1:
        return 1
    
    p_factors_count = 0
    temp_m = m
    i = 2
    while i * i <= temp_m:
        if temp_m % i == 0:
            p_factors_count += 1
            temp_m //= i
            # If the same prime factor divides again, m is not square-free
            if temp_m % i == 0:
                return 0
        i += 1
    if temp_m > 1:  # The remaining number is a prime factor
        p_factors_count += 1
        
    return -1 if p_factors_count % 2 != 0 else 1

def get_divisors(n):
    """Returns a sorted list of all positive divisors of n."""
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n // i)
    return sorted(list(divs))

def solve_bch_count():
    """
    Calculates the number of nonzero BCH coefficients for a given order n.
    """
    n = 10  # Order of the expansion
    k = 2   # Number of generators (X and Y)

    divisors = get_divisors(n)
    
    total_sum = 0
    calculation_terms = []

    # Apply Witt's formula: (1/n) * sum_{d|n} mu(n/d) * k^d
    for d in divisors:
        mu_val = moebius(n // d)
        term = mu_val * (k**d)
        total_sum += term
        calculation_terms.append(str(term))
    
    result = total_sum // n
    
    # Format the equation string for printing
    # Example: (1/10) * (2 - 4 - 32 + 1024) = 99
    sum_string = " + ".join(calculation_terms).replace("+ -", "- ")
    
    print(f"The number of terms is given by Witt's formula L_n(k) = (1/n) * Sum_{{d|n}} mu(n/d) * k^d")
    print(f"For n=10, k=2:")
    print(f"The divisors of 10 are: {divisors}")
    print("\nCalculating the sum:")
    sum_calc_details = []
    for d in divisors:
        sum_calc_details.append(f"d={d}: mu(10/{d})*2^{d} = mu({10//d})*2^{d} = {moebius(10//d)}*{2**d} = {moebius(10//d) * (2**d)}")
    print("\n".join(sum_calc_details))
    
    print("\nFinal Equation:")
    print(f"(1/{n}) * ({sum_string}) = (1/{n}) * ({total_sum}) = {result}")

if __name__ == "__main__":
    solve_bch_count()
    print("\n<<<99>>>")
