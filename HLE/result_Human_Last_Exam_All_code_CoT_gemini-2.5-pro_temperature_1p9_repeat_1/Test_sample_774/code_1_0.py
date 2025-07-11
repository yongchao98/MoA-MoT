import sympy
import numpy as np

def solve_problem():
    """
    Solves the given mathematical problem step-by-step.
    """

    # Part 1: The Summation
    # The term l(n,p) is the injectivity radius of the Stiefel manifold with the
    # Euclidean metric. This is known to be pi, independent of n and p.
    # The summation is over a 10x10 grid.
    sum_val = 10 * 10 * np.pi
    
    # Part 2: The Integral
    # The integral splits into two parts.
    # The second part simplifies to integral of x*exp(-x) from 0 to infinity, which is 1.
    integral_part_2 = 1
    
    # The first part depends on the dimensions of two manifolds.
    # Let's compute these dimensions.
    # Dimension formula: dim(M(n,p)) = n*p - p*(p+1)/2
    
    # Primes for d1
    n1_idx = 8231
    p1_idx = 781
    n1 = sympy.prime(n1_idx)
    p1 = sympy.prime(p1_idx)
    
    # Primes for d2
    n2_idx = 10231
    p2_idx = 2321
    n2 = sympy.prime(n2_idx)
    p2 = sympy.prime(p2_idx)
    
    # Calculate dimensions d1 and d2
    d1 = n1 * p1 - (p1 * (p1 + 1)) // 2
    d2 = n2 * p2 - (p2 * (p2 + 1)) // 2
    
    # Because d1 and d2 are extremely large, the first part of the integral
    # evaluates to 0. The integrand itself tends to 0 everywhere.
    integral_part_1 = 0
    
    total_integral = integral_part_1 + integral_part_2
    
    # Final Result
    final_result = sum_val * total_integral
    
    # Print the equation with all the computed numbers
    print("This problem computes the expression (Sum) * (Integral).")
    print("\n--- Part 1: Sum ---")
    print(f"l(n,p) is the injectivity radius, which equals pi.")
    print(f"The sum is 10 * 10 * pi = {sum_val:.10f}")
    
    print("\n--- Part 2: Integral ---")
    print(f"The integral is composed of two parts, I_1 + I_2.")
    print(f"I_2 = integral(x*exp(-x) dx) from 0 to inf = {integral_part_2}")
    
    print("\nI_1 depends on manifold dimensions d1 and d2:")
    print(f"n1 = p_({n1_idx}) = {n1}")
    print(f"p1 = p_({p1_idx}) = {p1}")
    print(f"d1 = {n1} * {p1} - ({p1} * ({p1} + 1)) / 2 = {d1}")
    
    print(f"\nn2 = p_({n2_idx}) = {n2}")
    print(f"p2 = p_({p2_idx}) = {p2}")
    print(f"d2 = {n2} * {p2} - ({p2} * ({p2} + 1)) / 2 = {d2}")
    
    print(f"\nSince d1 and d2 are very large, the integral I_1 evaluates to {integral_part_1}.")
    print(f"The total integral is I_1 + I_2 = {integral_part_1} + {integral_part_2} = {total_integral}")
    
    print("\n--- Final Calculation ---")
    print(f"Result = (Sum) * (Integral)")
    print(f"Result = ({sum_val:.10f}) * ({total_integral})")
    print(f"Final Value = {final_result:.10f}")

solve_problem()