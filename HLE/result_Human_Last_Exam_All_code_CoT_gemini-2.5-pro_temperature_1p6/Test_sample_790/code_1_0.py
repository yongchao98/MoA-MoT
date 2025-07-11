import math

class TitanFraction:
    """
    A class to represent and operate on fractions under Titan's 5-bit constraints.
    Note: In a real Titan simulation, checks would be added to ensure n and d
    are always <= 31. For this script, we ensure this through our calculation logic.
    """
    def __init__(self, numerator, denominator):
        self.n = numerator
        self.d = denominator

    def __str__(self):
        return f"{self.n}/{self.d}"

    def value(self):
        return self.n / self.d

def titan_multiply(frac1, frac2):
    """
    Multiplies two TitanFractions, simplifying before multiplication
    to adhere to the 5-bit constraint on intermediate values.
    This function finds the greatest common divisor (GCD) to simplify.
    """
    # Simplify frac1.n and frac2.d
    common1 = math.gcd(frac1.n, frac2.d)
    n1 = frac1.n // common1
    d2 = frac2.d // common1

    # Simplify frac2.n and frac1.d
    common2 = math.gcd(frac2.n, frac1.d)
    n2 = frac2.n // common2
    d1 = frac1.d // common2

    # The new numerators and denominators are pre-simplified.
    # The multiplication n1*n2 and d1*d2 is the final operation.
    final_n = n1 * n2
    final_d = d1 * d2

    return TitanFraction(final_n, final_d)

def solve_monkey_problem():
    """
    Solves the problem using Titan's computational rules.
    """
    print("Monkey and Coconut Problem on the Titan Supercomputer")
    print("----------------------------------------------------\n")
    print("The derived formula for the required force is: F = 0.3 * pi * g * sqrt(2)\n")

    # Step 1: Choose fractional approximations for all constants.
    # These values are chosen to allow for simplification and minimize final error.
    c_0_3 = TitanFraction(3, 10)
    # Using pi ~= 10/3 (3.33...) has a larger error but allows for perfect cancellation.
    pi_approx = TitanFraction(10, 3)
    # Using g = 10 is a common and simple approximation.
    g_approx = TitanFraction(10, 1)
    # Using sqrt(2) ~= 13/10 (1.3) allows for perfect cancellation.
    sqrt2_approx = TitanFraction(13, 10)

    print("Chosen approximations to satisfy 5-bit (0-31) constraints:")
    print(f"  0.3      ≈ {c_0_3}")
    print(f"  pi       ≈ {pi_approx}")
    print(f"  g        ≈ {g_approx}")
    print(f"  sqrt(2)  ≈ {sqrt2_approx}\n")

    # Step 2: Calculate F step-by-step
    print("Calculation Steps:")
    
    # Calculate term1 = 0.3 * pi
    print("1. Calculate first part: (0.3 * pi)")
    print(f"   {c_0_3} * {pi_approx}")
    term1 = titan_multiply(c_0_3, pi_approx)
    print(f"   = {term1.n}/{term1.d}  (simplified from 30/30)")
    
    # Calculate term2 = g * sqrt(2)
    print("\n2. Calculate second part: (g * sqrt(2))")
    print(f"   {g_approx} * {sqrt2_approx}")
    term2 = titan_multiply(g_approx, sqrt2_approx)
    print(f"   = {term2.n}/{term2.d}  (simplified from 130/10 by cancelling the 10s before multiplying)")

    # Final multiplication
    print("\n3. Calculate the final force F by multiplying the two parts:")
    print(f"   F = {term1} * {term2}")
    final_force_frac = titan_multiply(term1, term2)
    print(f"   = {final_force_frac}")
    
    # Final Result and Error Analysis
    final_force_val = final_force_frac.value()
    print(f"\nThe final calculated force is {final_force_frac}, which is {final_force_val:.1f} Newtons.")

    # Calculate error against the more precise theoretical value
    true_force_val = 0.3 * math.pi * 9.8 * math.sqrt(2)
    error = abs(final_force_val - true_force_val)
    
    print(f"\nThe theoretical force is ~{true_force_val:.4f} N.")
    print(f"The absolute error of our calculation is {error:.3f}.\n")
    
    # Check if the calculated force hits the target (9.9m to 10.1m height)
    # F_target = F_center * (0.5 / (1 - y/x))
    # Target range for y is 9.9 to 10.1
    f_low = true_force_val * (0.5 / (1 - 9.9/20))
    f_high = true_force_val * (0.5 / (1 - 10.1/20))
    
    print(f"For the rock to hit the coconut, the force must be between ~{f_low:.2f} N and ~{f_high:.2f} N.")
    
    if f_low <= final_force_val <= f_high:
        print(f"Our calculated force of {final_force_val:.1f} N is within the required range. The calculation is successful.")
        # Return the final answer in the requested format
        print(f"\n<<<Y[{error:.3f}]>>>")
    else:
        print(f"Our calculated force of {final_force_val:.1f} N is NOT within the required range. The calculation failed.")
        print("\n<<<N0>>>")


if __name__ == '__main__':
    solve_monkey_problem()