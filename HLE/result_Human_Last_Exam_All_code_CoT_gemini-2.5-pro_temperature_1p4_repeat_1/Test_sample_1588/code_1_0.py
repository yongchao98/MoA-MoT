import math

def get_prime_factorization(num):
    """Returns a dictionary of prime factors and their powers."""
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

def mobius(num):
    """Calculates the Mobius function mu(num)."""
    if num == 1:
        return 1
    
    factors = get_prime_factorization(num)
    for p in factors:
        if factors[p] > 1:
            return 0
            
    return (-1) ** len(factors)

def solve_bch_order_10():
    """
    Calculates the number of nonzero coefficients of order 10 in the BCH expansion
    using Witt's formula.
    """
    n = 10
    k = 2

    # Find divisors of n
    divisors = sorted([d for d in range(1, n + 1) if n % d == 0])

    numerator_sum = 0
    equation_parts = []
    calculation_parts = []

    for d in divisors:
        mu_d = mobius(d)
        power = n // d
        term_val = mu_d * (k ** power)
        numerator_sum += term_val
        
        # Build strings for printing the equation and calculation
        part_str = f"μ({d}) * {k}^({n}/{d})"
        equation_parts.append(part_str)
        
        calc_str = f"({mu_d}) * {k}^{power}"
        calculation_parts.append(calc_str)

    result = numerator_sum // n

    # Print the explanation and final equation
    print("The number of nonzero coefficients is given by Witt's formula for the dimension of the free Lie algebra L_k(n):")
    print(f"L_k(n) = (1/n) * Σ_{{d|n}} μ(d) * k^(n/d)")
    print("\nFor this problem, n=10 and k=2. The divisors of 10 are 1, 2, 5, 10.")
    
    print("\nStep 1: Write out the formula for n=10, k=2")
    final_equation = f"L_2(10) = (1/{n}) * [ " + " + ".join(equation_parts) + " ]"
    print(final_equation)
    
    print("\nStep 2: Substitute the values of the Mobius function and powers")
    calculation_str = f"L_2(10) = (1/{n}) * [ " + " + ".join(calculation_parts) + " ]"
    print(calculation_str)

    print("\nStep 3: Calculate the terms inside the brackets")
    final_calculation = f"L_2(10) = (1/{n}) * [ {numerator_sum} ]"
    print(final_calculation)
    
    print("\nStep 4: Final result")
    print(f"L_2(10) = {result}")


solve_bch_order_10()
<<<99>>>