import math
from fractions import Fraction

def solve_problem():
    """
    Solves the physics problem based on Wuxing architecture constraints.
    """
    print("1. Deriving the governing equation from physics principles.")
    print("   The collision condition leads to a quadratic equation for the initial speed u:")
    print("   (sin(2*a)) * u^2 + (2*v*sin(a)) * u - (d*g) = 0")
    print("\n2. Adapting for Wuxing architecture.")
    print("   - Using fractional approximations for constants.")
    print("   - Using Newton-Raphson method to solve the equation without sqrt().")
    
    # Constants for the problem
    v = 5       # lion's speed in m/s
    d = 300     # initial distance in m
    a_deg = 60  # throw angle in degrees

    # Wuxing-compatible fractional approximations for physical/math constants
    # g = 9.8 = 98/10
    g = Fraction(98, 10)
    # sin(60) = sin(120) approx 0.86 = 86/100 = 43/50
    # This approximation avoids overflow on the Wuxing architecture's 'frac' type.
    sin_a = Fraction(43, 50)
    sin_2a = sin_a # Since sin(120) is also approx 0.86
    
    # The quadratic equation is A*u^2 + B*u + C = 0
    # Calculate coefficients A, B, C using fraction arithmetic
    A = sin_2a
    B = 2 * v * sin_a
    C = -d * g

    # As required, print each number in the final equation
    print("\n3. The final equation to be solved:")
    print(f"   ({A.numerator}/{A.denominator}) * u^2 + ({B.numerator}/{B.denominator}) * u + ({C.numerator}) = 0")
    
    # 4. Solve for u using Newton-Raphson method
    # f(u) = A*u^2 + B*u + C
    # f'(u) = 2*A*u + B
    u = Fraction(50)  # Initial guess
    
    # Iteratively find the root
    for _ in range(5):  # 5 iterations are sufficient for convergence
        f_u = A * u**2 + B * u + C
        f_prime_u = 2 * A * u + B
        if f_prime_u == 0:
            break
        u = u - f_u / f_prime_u
        
    # Round final result to one decimal place
    u_final = round(float(u), 1)

    # 5. Calculate memory usage 'z' for the variables
    # frac A (6D) + frac B (6D) + int C (5D) + frac u (6D)
    # A, B, C, and u are the main variables stored.
    mem_A = 6  # frac
    mem_B = 6  # frac
    mem_C = 5  # int
    mem_u = 6  # frac
    z = mem_A + mem_B + mem_C + mem_u
    print(f"\n4. Calculating memory usage (z):")
    print(f"   Memory for A(frac), B(frac), C(int), u(frac) is {mem_A}+{mem_B}+{mem_C}+{mem_u} = {z} D.")

    # 6. Final answer formatting
    result = f"{u_final}:{z}"
    print("\n---")
    print(f"Final calculated speed u = {u_final} m/s")
    print(f"Final calculated memory z = {z} D")
    print(f"\n<<<Formatted Answer>>>")
    print(f"<<<{result}>>>")

solve_problem()