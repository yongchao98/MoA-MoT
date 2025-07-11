import math

def get_prime_factorization(num):
    """
    Returns a dictionary of prime factors of num.
    e.g., get_prime_factorization(12) -> {2: 2, 3: 1}
    """
    factors = {}
    d = 2
    temp = num
    while d * d <= temp:
        while temp % d == 0:
            factors[d] = factors.get(d, 0) + 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] = factors.get(temp, 0) + 1
    return factors

def mobius(n):
    """
    Calculates the Mobius function mu(n).
    mu(n) = 1 if n is a square-free with an even number of prime factors.
    mu(n) = -1 if n is a square-free with an odd number of prime factors.
    mu(n) = 0 if n has a squared prime factor.
    """
    if not isinstance(n, int) or n < 1:
        raise ValueError("Input must be a positive integer.")
    if n == 1:
        return 1
    
    factors = get_prime_factorization(n)
    for p in factors:
        if factors[p] > 1:
            return 0 # Has a squared prime factor
            
    return (-1)**len(factors) # (-1)^k for k distinct prime factors

def get_divisors(n):
    """
    Returns a sorted list of all positive divisors of n.
    """
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n//i)
    return sorted(list(divs))

def solve_bch_nonzero_coefficients():
    """
    Calculates the number of nonzero coefficients of a given order in the
    BCH expansion using Witt's formula and prints the calculation steps.
    """
    # 1. Define parameters for the problem
    n = 10  # The order of the coefficients
    k = 2   # The number of generators (X and Y)

    print("To find the number of nonzero coefficients of order 10 in the Baker-Campbell-Hausdorff expansion,")
    print("we use Witt's formula for the dimension of the degree-n subspace of a free Lie algebra on k generators.\n")
    
    print("Witt's Formula: L_k(n) = (1/n) * sum_{d|n} mu(d) * k^(n/d)")
    print(f"In this case, n = {n} and k = {k}.\n")

    # 2. Get divisors of n
    divisors = get_divisors(n)
    print(f"The positive divisors 'd' of {n} are: {divisors}\n")

    # 3. Calculate each term in the sum
    total_sum = 0
    calculation_parts = []
    
    print("Calculating each term mu(d) * k^(n/d) for each divisor d:")
    for d in divisors:
        mu_d = mobius(d)
        power_val = n // d
        term_result = mu_d * (k ** power_val)
        
        print(f"For d={d}: mu({d}) * {k}^({n}/{d}) = {mu_d} * {k}^{power_val} = {term_result}")
        
        total_sum += term_result
        calculation_parts.append(term_result)

    print("\nNow, we construct the full equation:")
    
    # 4. Assemble the final equation string for printing
    # Create the string for the sum, e.g., "1024 - 32 - 4 + 2"
    sum_string = ""
    for i, term in enumerate(calculation_parts):
        if i == 0:
            sum_string += f"{term}"
        elif term >= 0:
            sum_string += f" + {term}"
        else:
            sum_string += f" - {abs(term)}"
            
    print(f"L_{k}({n}) = (1/{n}) * ( {sum_string} )")
    print(f"L_{k}({n}) = (1/{n}) * ( {total_sum} )")
    
    # 5. Calculate final result
    final_result = total_sum // n
    
    print(f"L_{k}({n}) = {final_result}\n")
    print("----------------------------------------------------------------")
    print(f"The number of nonzero coefficients of order 10 is {final_result}.")
    print("----------------------------------------------------------------")

if __name__ == "__main__":
    solve_bch_nonzero_coefficients()