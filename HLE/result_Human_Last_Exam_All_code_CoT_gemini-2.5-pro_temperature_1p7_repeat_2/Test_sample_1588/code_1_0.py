import math

def get_prime_factorization(num):
    """Returns a dictionary of prime factors and their powers."""
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
    """Calculates the Mobius function mu(n)."""
    if n == 1:
        return 1
    factors = get_prime_factorization(n)
    for p in factors:
        if factors[p] > 1:  # n has a squared prime factor
            return 0
    if len(factors) % 2 != 0:  # n has an odd number of distinct prime factors
        return -1
    else:  # n has an even number of distinct prime factors
        return 1

def main():
    """
    Calculates the number of nonzero BCH coefficients of order 10
    using Witt's formula.
    """
    n = 10
    k = 2

    # Find all positive divisors of n
    divisors = [d for d in range(1, n + 1) if n % d == 0]

    sum_terms_values = []
    sum_terms_str_parts = []
    
    # Calculate each term of the sum in Witt's formula
    for d in divisors:
        mu_d = mobius(d)
        power = n // d
        term_val = mu_d * (k ** power)
        
        sum_terms_values.append(term_val)
        sum_terms_str_parts.append(f"({mu_d}) * {k}^{power}")

    # Calculate the final result
    total_sum = sum(sum_terms_values)
    result = total_sum / n
    
    print("The number of nonzero coefficients is calculated using Witt's formula:")
    print(f"d(n) = (1/n) * sum_{{d|n}} (mu(d) * k^(n/d))")
    print("\nFor n=10, k=2:")
    
    # Print the equation with all numbers filled in
    equation_str = f"d(10) = (1/{n}) * [ " + " + ".join(sum_terms_str_parts) + " ]"
    print(equation_str)
    
    # Print the equation with evaluated terms
    evaluated_terms = " + ".join(map(str, sum_terms_values)).replace('+ -', '- ')
    evaluated_str = f"d(10) = (1/{n}) * [ {evaluated_terms} ]"
    print(evaluated_str)
    
    # Print the final calculation step and result
    final_calc_str = f"d(10) = (1/{n}) * {total_sum} = {int(result)}"
    print(final_calc_str)
    
    print(f"\nFinal Answer: {int(result)}")

if __name__ == "__main__":
    main()