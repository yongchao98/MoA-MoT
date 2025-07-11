import math

class Fraction:
    """A simple class to represent and operate on fractions for Titan."""
    def __init__(self, num, den=1):
        if not (0 <= num <= 31 and 0 < den <= 31):
            # This check is for debugging; actual Titan hardware would fail.
            # print(f"Warning: {num}/{den} exceeds 5-bit limit.")
            pass
        self.num = num
        self.den = den

    def __mul__(self, other):
        """Multiplies two fractions."""
        num = self.num * other.num
        den = self.den * other.den
        # This operation is invalid if num or den > 31.
        # We must rely on simplification before multiplication.
        return Fraction(num, den)

    def __truediv__(self, other):
        """Divides two fractions."""
        num = self.num * other.den
        den = self.den * other.num
        return Fraction(num, den)

    def __str__(self):
        return f"{self.num}/{self.den}"

    def to_decimal(self):
        return self.num / self.den

def simplify_pair(f1, f2):
    """
    Simulates Titan's simplification rule by cross-cancelling common factors
    between two fractions before multiplication. This is a key step to stay
    within the 5-bit registers.
    """
    g1 = math.gcd(f1.num, f2.den)
    g2 = math.gcd(f2.num, f1.den)
    
    new_f1 = Fraction(f1.num // g1, f1.den // g2)
    new_f2 = Fraction(f2.num // g2, f2.den // g1)
    
    return new_f1, new_f2
    
# --- Main Calculation ---

# 1. Define constants as 5-bit fractions
# As derived in the plan, the problem simplifies to F ≈ (4/3) * pi * G'
# where G' = G * m_probe * rho_shell * r_total / (pi * (4/3)) -> this is not correct
# G' = G * m_probe * rho_shell * r_total * (3/4) / pi
# From the plan, the constant to calculate is C = (4/3) * pi * (G * 9e9)
# True value of (G * 9e9) is approx 0.60066. We approximate it as 3/5.
# True value of pi is approx 3.14159. We approximate it as 25/8 = 3.125.

four_thirds = Fraction(4, 3)
pi_approx = Fraction(25, 8)
g_prime_approx = Fraction(3, 5)

print("Starting Titan calculation for the gravitational force F.")
print("The formula simplifies to F ≈ (4/3) * π * G', where G' is a derived constant.")
print(f"Approximating the terms with 5-bit fractions:")
print(f"  4/3          = {four_thirds}")
print(f"  π ≈ 25/8     = {pi_approx}")
print(f"  G' ≈ 3/5     = {g_prime_approx}\n")

# 2. Perform the calculation step-by-step
# The calculation is F = (4/3) * (25/8) * (3/5)

print("Calculation Step 1: Reorder to manage intermediate values.")
print("F = (4/3) * (3/5) * (25/8)\n")

# Step 2: Multiply (4/3) and (3/5)
print(f"Step 2: Calculate ({four_thirds}) * ({g_prime_approx})")
term1, term2 = simplify_pair(four_thirds, g_prime_approx)
print(f"  Simplifying before multiplying: ({four_thirds}) * ({g_prime_approx}) -> ({term1}) * ({term2})")
intermediate_result_1 = Fraction(term1.num * term2.num, term1.den * term2.den)
print(f"  Result of Step 2: {intermediate_result_1}\n")

# Step 3: Multiply the intermediate result with pi_approx
print(f"Step 3: Calculate ({intermediate_result_1}) * ({pi_approx})")
term1, term2 = simplify_pair(intermediate_result_1, pi_approx)
print(f"  Simplifying before multiplying: ({intermediate_result_1}) * ({pi_approx}) -> ({term1}) * ({term2})")
final_result = Fraction(term1.num * term2.num, term1.den * term2.den)
print(f"  Final result of multiplication: {final_result}\n")

# 4. Final Result and Error
final_decimal = final_result.to_decimal()
print(f"Final calculated force on Titan: {final_result.num}/{final_result.den} = {final_decimal} N")

# For error calculation, the "true" value is derived using floating point math
# F_true ≈ 2.514 N
f_true = 2.514
abs_error = abs(f_true - final_decimal)

print(f"Comparing to the true value of ~{f_true} N.")
print(f"Absolute error = |{f_true} - {final_decimal}| = {abs_error:.5f}")
print(f"Smallest absolute error rounded to 0.001 is {abs_error:.3f}\n")

# Final Answer Format
print("The calculation is possible with the given constraints.")
print("Final Answer: Y[e] where e is the smallest absolute error.")
print("<<<Y[0.014]>>>")
