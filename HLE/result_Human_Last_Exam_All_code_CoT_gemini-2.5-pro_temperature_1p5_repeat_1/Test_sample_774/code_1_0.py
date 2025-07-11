import math
import sympy

def solve():
    """
    Solves the given mathematical problem by breaking it down into a sum and an integral,
    analyzing each part, and then combining the results.
    """
    
    # Part 1: The Summation
    # The term l(n,p) represents the injectivity radius of the Stiefel manifold M(n,p)
    # with the standard metric. This is a known result from differential geometry, l(n,p) = pi.
    # The problem specifies n, p >= 5. We must verify that for the given prime indices,
    # n >= p.
    # The sum is Sum_{i=1 to 10} Sum_{j=1 to 10} l(p_{21367+i}, p_{14567+j}).
    # Smallest n is p_{21367+1} = p_{21368}.
    # Largest p is p_{14567+10} = p_{14577}.
    # Since 21368 > 14577, n > p is always satisfied.
    # All primes with these indices are much larger than 5.
    # Therefore, every term in the sum is pi.
    
    num_i = 10
    num_j = 10
    l_val = math.pi
    sum_val = num_i * num_j * l_val
    
    print("--- Part 1: Calculation of the Sum ---")
    print(f"The value of the injectivity radius l(n,p) is a constant: {l_val}")
    print(f"The sum is over i from 1 to {num_i} and j from 1 to {num_j}.")
    print(f"The total sum is {num_i} * {num_j} * {l_val} = {sum_val}")
    print("-" * 35)
    print()

    # Part 2: The Integral
    # The integral is simplified by splitting the integrand.
    # Integral = Integral[ (term1) + (term2) ]dx
    # term1 = (x^(2d1) - x^(2d2)) / (x * (1+x^(2d1)) * (1+x^(2d2)) * sqrt(e^(2x)-1))
    # term2 = x * exp(-x)
    # The integral of term1 approaches 0 because the dimensions d1 and d2 are very large.
    # The integral of term2 from 0 to infinity is exactly 1.
    # First, we calculate d1 and d2 to confirm they are large.
    
    # dim(M(n,p)) = n*p - p*(p+1)/2
    
    n1_idx = 8231
    p1_idx = 781
    n1 = sympy.prime(n1_idx)
    p1 = sympy.prime(p1_idx)
    d1 = n1 * p1 - (p1 * (p1 + 1)) // 2

    n2_idx = 10231
    p2_idx = 2321
    n2 = sympy.prime(n2_idx)
    p2 = sympy.prime(p2_idx)
    d2 = n2 * p2 - (p2 * (p2 + 1)) // 2

    # Because d1 and d2 are large, the first part of the integral tends to 0.
    integral_part1 = 0
    # The second part of the integral, Integral[x*e^(-x) dx] from 0 to inf, is 1.
    integral_part2 = 1
    integral_val = integral_part1 + integral_part2

    print("--- Part 2: Calculation of the Integral ---")
    print("Numbers for dimension d1:")
    print(f"n1 = p_({n1_idx}) = {n1}")
    print(f"p1 = p_({p1_idx}) = {p1}")
    print(f"Dimension d1 = {n1}*{p1} - {p1}*({p1}+1)/2 = {d1}")
    print("\nNumbers for dimension d2:")
    print(f"n2 = p_({n2_idx}) = {n2}")
    print(f"p2 = p_({p2_idx}) = {p2}")
    print(f"Dimension d2 = {n2}*{p2} - {p2}*({p2}+1)/2 = {d2}")
    
    print(f"\nAs d1 and d2 are very large, the first component of the integral evaluates to {integral_part1}.")
    print(f"The second component of the integral evaluates to {integral_part2}.")
    print(f"The total value of the integral is {integral_part1} + {integral_part2} = {integral_val}")
    print("-" * 35)
    print()

    # Final Result
    final_result = sum_val * integral_val
    
    print("--- Final Calculation ---")
    print("The final equation is: (Sum Value) * (Integral Value)")
    print(f"Final Value = {sum_val} * {integral_val}")
    print(f"The final calculated value is: {final_result}")

solve()