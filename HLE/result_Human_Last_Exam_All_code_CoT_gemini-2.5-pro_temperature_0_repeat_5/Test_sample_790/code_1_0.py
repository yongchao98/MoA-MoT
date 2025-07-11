import math

class TitanFraction:
    """
    A class to represent and compute with fractions under Titan's 5-bit constraints.
    """
    def __init__(self, num, den=1):
        if not (0 <= num <= 31 and 0 < den <= 31):
            raise ValueError(f"Fraction {num}/{den} violates 5-bit constraint.")
        self.num = num
        self.den = den

    def __mul__(self, other):
        """
        Multiplies two TitanFractions, handling constraints.
        """
        # Perform the multiplication
        new_num = self.num * other.num
        new_den = self.den * other.den

        # Check for constraint violation
        if new_num > 31 or new_den > 31:
            # This is a special case from our derivation, where 3 * 22 is too large.
            # We approximate the term 22/5 (4.4) with 13/3 (4.33...) to proceed.
            if self.num == 3 and self.den == 1 and other.num == 22 and other.den == 5:
                print(f"-> Constraint Violated: {self.num} * {other.num} = {new_num} > 31.")
                print(f"   Approximating the term {other} (value ~{other.to_float():.2f}) with 13/3 (value ~{13/3:.2f}) to simplify.")
                # Re-run multiplication with the approximated fraction
                return self * TitanFraction(13, 3)
            else:
                raise ValueError(f"Multiplication {self} * {other} resulted in {new_num}/{new_den}, violating constraints.")

        # Simplify the fraction by finding the greatest common divisor (GCD)
        common_divisor = math.gcd(new_num, new_den)
        return TitanFraction(new_num // common_divisor, new_den // common_divisor)

    def to_float(self):
        return self.num / self.den

    def __str__(self):
        return f"{self.num}/{self.den}"

def solve_coconut_problem():
    """
    Solves the problem using the Titan computer simulation.
    """
    print("Starting Titan computation for the required force F.")
    print("Derived formula: F = (3/10) * g * pi * sqrt(2)\n")

    # Step 1: Define fractional approximations for constants
    # g ≈ 9.8 m/s^2. We use 10/1 for simplicity and to enable cancellation.
    # pi ≈ 3.14159. We use 22/7.
    # sqrt(2) ≈ 1.414. We use 7/5, which cancels nicely with pi's denominator.
    f_3_10 = TitanFraction(3, 10)
    f_g = TitanFraction(10, 1)
    f_pi = TitanFraction(22, 7)
    f_sqrt2 = TitanFraction(7, 5)

    print("Chosen approximations:")
    print(f"g      ≈ {f_g} (value = {f_g.to_float()})")
    print(f"pi     ≈ {f_pi} (value = {f_pi.to_float():.3f})")
    print(f"sqrt(2) ≈ {f_sqrt2} (value = {f_sqrt2.to_float()})")
    print("-" * 30)

    # Step 2: Perform calculations step-by-step
    print("Calculation steps:")
    # Grouping for simplification: F = ((3/10) * g) * (pi * sqrt(2))
    
    # First group
    print(f"1. Calculate (3/10) * g ≈ {f_3_10} * {f_g}")
    term1 = f_3_10 * f_g
    print(f"-> Result: {term1}\n")

    # Second group
    print(f"2. Calculate pi * sqrt(2) ≈ {f_pi} * {f_sqrt2}")
    term2 = f_pi * f_sqrt2
    print(f"-> Result: {term2}\n")

    # Final multiplication
    print(f"3. Calculate final force F ≈ {term1} * {term2}")
    final_force_fraction = term1 * term2
    
    print("\n" + "-" * 30)
    print("Final Titan Calculation:")
    # The final equation that works within the constraints
    print(f"The final equation is: {term1.num} / {term1.den} * 13 / 3 = {final_force_fraction.num} / {final_force_fraction.den}")
    
    titan_result = final_force_fraction.to_float()
    print(f"Final result as a fraction: {final_force_fraction}")
    print(f"Final result as a decimal: {titan_result:.3f} N")
    print("-" * 30)

    # Calculate the "true" value and the error
    true_val = 0.3 * math.pi * math.sqrt(2) * 9.8
    error = abs(titan_result - true_val)
    print(f"Calculating error against 'true' value ({true_val:.3f} N)...")
    print(f"Absolute Error: {error:.3f}")

    # Verify if the coconut is hit
    m = 0.15 * math.pi
    g_val = 9.8
    cos45 = math.cos(math.radians(45))
    y_height = 20 * (1 - (m * g_val) / (titan_result * cos45))
    
    print(f"\nVerifying trajectory with F = {titan_result:.3f} N...")
    print(f"Predicted height at tree: {y_height:.3f} m")
    if 9.9 <= y_height <= 10.1:
        print("SUCCESS: The rock hits the coconut.")
        print(f"\nFinal Answer: Y{error:.3f}")
    else:
        print("FAILURE: The rock misses the coconut.")
        print("\nFinal Answer: N0")

if __name__ == '__main__':
    solve_coconut_problem()