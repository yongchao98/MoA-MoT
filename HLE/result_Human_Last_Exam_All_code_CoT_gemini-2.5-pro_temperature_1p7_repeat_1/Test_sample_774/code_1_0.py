import math
import sympy

def solve_problem():
    """
    Solves the given mathematical problem by calculating the sum and integral parts.
    """
    
    # Part 1: Calculate the Sum
    # The term l(n,p) is the injectivity radius of the Stiefel manifold M(n,p),
    # which is known to be pi.
    # The sum is over i=1..10 and j=1..10.
    sum_val = 100 * math.pi
    
    # Part 2: Calculate the Integral
    # The integral splits into two parts, I = I_1 + I_2.
    # I_2 = integral from 0 to inf of x*exp(-x) dx = 1.
    
    # For I_1, we need the dimensions d1 and d2.
    # First, find the required prime numbers.
    # The k-th prime is given by sympy.prime(k)
    p_781 = sympy.prime(781)
    p_8231 = sympy.prime(8231)
    p_10231 = sympy.prime(10231)
    p_2321 = sympy.prime(2321)
    
    # Calculate dimension d1 = dim(M(p_8231, p_781))
    n1 = p_8231
    p1 = p_781
    d1 = n1 * p1 - (p1 * (p1 + 1)) // 2
    
    # Calculate dimension d2 = dim(M(p_10231, p_2321))
    n2 = p_10231
    p2 = p_2321
    d2 = n2 * p2 - (p2 * (p2 + 1)) // 2
    
    print(f"The dimension d1 is calculated from n=p_8231={n1} and p=p_781={p1}.")
    print(f"d1 = {n1} * {p1} - {p1} * ({p1} + 1) / 2 = {d1}")
    print(f"The dimension d2 is calculated from n=p_10231={n2} and p=p_2321={p2}.")
    print(f"d2 = {n2} * {p2} - {p2} * ({p2} + 1) / 2 = {d2}")

    # The integral I_1 has an integrand that tends to 0 everywhere because d1 and d2
    # are very large. Thus, I_1 is 0.
    I1 = 0
    I2 = 1  # integral of x*exp(-x) from 0 to inf
    integral_val = I1 + I2
    
    # Final Result
    final_result = sum_val * integral_val

    print("\n--- Final Equation Numbers ---")
    print(f"Value of the sum part: 100 * pi = {sum_val}")
    print(f"Value of the integral part: {integral_val}")
    print(f"Final Result = (Sum) * (Integral) = {sum_val} * {integral_val} = {final_result}")

solve_problem()