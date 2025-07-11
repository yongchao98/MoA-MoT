import math

def get_prime_factorization(num):
    """
    Returns a dictionary of prime factors of num and their powers.
    e.g., get_prime_factorization(12) returns {2: 2, 3: 1}
    """
    factors = {}
    d = 2
    temp = num
    while d * d <= temp:
        while (temp % d) == 0:
            factors[d] = factors.get(d, 0) + 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] = factors.get(temp, 0) + 1
    return factors

def mobius(n):
    """
    Calculates the Mobius function mu(n).
    mu(1) = 1
    mu(n) = 0 if n has a squared prime factor.
    mu(n) = (-1)^k if n is a product of k distinct primes.
    """
    if not isinstance(n, int) or n < 1:
        raise ValueError("Möbius function is defined only for positive integers.")
    if n == 1:
        return 1
    
    factors = get_prime_factorization(n)
    for p in factors:
        if factors[p] > 1:
            return 0
    return (-1)**len(factors)

def solve_bch_order_10():
    """
    Calculates the number of nonzero coefficients of order 10 in the BCH expansion.
    """
    n = 10
    k = 2
    
    # Find divisors of n
    divisors = sorted([i for i in range(1, n + 1) if n % i == 0])
    
    # Calculate each term in the Witt's formula sum
    total_sum = 0
    mu_values = []
    terms = []
    term_values = []
    
    for d in divisors:
        mu_val = mobius(n // d)
        mu_values.append(f"mu({n//d})={mu_val}")
        term_val = mu_val * (k**d)
        terms.append(f"({mu_val}) * {k}^{d}")
        term_values.append(str(term_val))
        total_sum += term_val
        
    result = total_sum // n
    
    # Print the explanation and the step-by-step calculation
    print(f"The number of nonzero coefficients is given by Witt's formula for n={n} and k={k}:")
    print(f"L_k(n) = (1/n) * sum_{{d|n}} mu(n/d) * k^d\n")
    
    print(f"For n={n}, the divisors 'd' are: {divisors}")
    print(f"The corresponding mu(n/d) values are: {', '.join(mu_values)}\n")
    
    print("Substituting these values into the formula:")
    
    equation_symbolic = f"L_2(10) = (1/{n}) * [ " + " + ".join(terms) + " ]"
    print(equation_symbolic.replace("+ -", "- "))

    equation_numeric = f"L_2(10) = (1/{n}) * [ " + " + ".join(term_values) + " ]"
    print(equation_numeric.replace("+ -", "- "))
    
    equation_summed = f"L_2(10) = (1/{n}) * [ {total_sum} ]"
    print(equation_summed)
    
    print(f"\nL_2(10) = {result}")
    
    print(f"\nThus, the number of nonzero coefficients of order 10 is {result}.")

if __name__ == "__main__":
    solve_bch_order_10()
