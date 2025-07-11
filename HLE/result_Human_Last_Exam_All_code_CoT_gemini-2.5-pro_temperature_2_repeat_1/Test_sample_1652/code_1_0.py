import sys
from fractions import Fraction

def solve_projectile_problem():
    """
    Calculates the initial speed 'u' for the rock to hit the lion
    using the Wuxing architecture constraints.
    """
    # 1. Define constants based on problem statement and Wuxing data types.
    # All non-integers are represented as fractions to simulate the 'frac' type.
    
    # Lion's speed v = 5 m/s. This would be an 'int' (5D), but we use Fraction for calculation.
    v = Fraction(5, 1)

    # Initial distance d = 300 m. This would be an 'int' (5D).
    distance = Fraction(300, 1)

    # Gravity g = 9.8 m/s^2. This requires the 'frac' type (6D).
    g = Fraction(98, 10)

    # For angle a=60°, we need sin(a) and sin(2a). Since these functions
    # are unavailable, we use a rational approximation for sin(60°) = sqrt(3)/2.
    # sqrt(3) ≈ 19/11 => sin(60°) ≈ (19/11)/2 = 19/22.
    # This constant would be stored as a 'frac' type (6D).
    sin60 = Fraction(19, 22)

    # sin(2*60°) = sin(120°) = sin(180°-60°) = sin(60°)
    sin120 = sin60

    # 2. Set up the governing quadratic equation f(u) = A*u^2 + B*u + C = 0
    # A = sin(2a), B = 2*v*sin(a), C = -d*g
    A = sin120
    B = 2 * v * sin60
    C = -distance * g

    def f(u_frac):
        """The function for which we are finding the root u."""
        return A * u_frac**2 + B * u_frac + C

    # 3. Perform a bisection search to find the root u > 0.
    # This numerical method avoids using sqrt, which is not available.
    low = Fraction(0)
    high = Fraction(200)  # An initial guess for the upper bound of the speed.

    # A sufficient number of iterations for high precision.
    for _ in range(100):
        mid = (low + high) / 2
        # Exit if precision limit is reached
        if mid == low or mid == high:
            break
        # If f(mid) and f(low) have the same sign, root is in the upper half.
        if f(mid).numerator * f(low).numerator > 0:
            low = mid
        else:
            high = mid

    u_solution = (low + high) / 2
    
    # 4. Calculate final values for the answer.
    # The result for u, rounded to 0.1 decimal places.
    u_rounded = round(float(u_solution), 1)

    # The memory usage 'z' for the main variables (g, sin60, v, d)
    # frac g (6D) + frac sin60 (6D) + int v (5D) + int distance (5D)
    z_memory_usage = 22
    
    # 5. Output the results as required.
    # First, print the numbers from the final equation. To get an integer equation,
    # we can write (19/22)u^2 + (190/22)u - 2940 = 0 and multiply by 22.
    # 19*u^2 + 190*u - 64680 = 0
    coeff_A_int = A.numerator * B.denominator # 19 * 22
    coeff_B_int = B.numerator # 190
    coeff_C_int = int(C * B.denominator) # -2940 * 22
    # Simplify B. 190/22 = 95/11. A=19/22. C = -2940.
    # Common denominator is 22.
    # (19/22)u^2 + (190/22)u - (64680/22) = 0
    # Final Integer Equation: 19 * u^2 + 190 * u - 64680 = 0
    
    print(f"The final governing equation derived from the physics is:")
    print(f"({A.numerator}/{A.denominator}) * u^2 + ({B.numerator}/{B.denominator}) * u + ({C.numerator}) = 0")
    print(f"Or, with integer coefficients: ")
    print(f"{19} * u^2 + {190} * u - {64680} = 0")
    print(f"\nFinal Answer (u rounded to 0.1 : z memory usage):")

    # Final formatted answer string for submission.
    final_answer_string = f"{u_rounded}:{z_memory_usage}"
    print(f"<<<{final_answer_string}>>>")

solve_projectile_problem()