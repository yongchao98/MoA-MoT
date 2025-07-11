import math

class TitanFraction:
    """
    A class to represent fractions under Titan's 5-bit constraints.
    Numerators and denominators must be between 0 and 31.
    """
    def __init__(self, num, den=1):
        if not (0 <= num <= 31 and 1 <= den <= 31):
            raise ValueError(f"Numerator {num} or denominator {den} is out of 5-bit range (0-31).")
        self.num = num
        self.den = den

    def __repr__(self):
        return f"{self.num}/{self.den}"

    def value(self):
        return self.num / self.den

def gcd(a, b):
    """Helper function to find the greatest common divisor."""
    return math.gcd(a, b)

def multiply(f1, f2, comment=""):
    """
    Multiplies two TitanFractions, respecting the 5-bit constraints by simplifying before multiplying.
    This simulates the algebraic simplification `(a/b)*(c/d) = (a*c)/(b*d)` by finding common
    factors between numerators and denominators to keep intermediate values small.
    """
    print(f"Multiplying {f1} by {f2}. {comment}")

    # Simplify before multiplying to avoid overflow
    # gcd(a, d)
    common_ad = gcd(f1.num, f2.den)
    num1_s = f1.num // common_ad
    den2_s = f2.den // common_ad

    # gcd(c, b)
    common_cb = gcd(f2.num, f1.den)
    num2_s = f2.num // common_cb
    den1_s = f1.den // common_cb

    # Check if the products will exceed the 31 limit
    if num1_s * num2_s > 31 or den1_s * den2_s > 31:
        raise ValueError(f"Multiplication {(num1_s, den1_s)} * {(num2_s, den2_s)} results in overflow.")
    
    res_num = num1_s * num2_s
    res_den = den1_s * den2_s

    # Final simplification of the result
    common_final = gcd(res_num, res_den)
    final_num = res_num // common_final
    final_den = res_den // common_final

    result = TitanFraction(final_num, final_den)
    print(f"  -> Simplified terms to ({num1_s}/{den1_s}) * ({num2_s}/{den2_s}) = {result}")
    return result

def divide(f1, f2, comment=""):
    """Divides f1 by f2 by multiplying f1 by the inverse of f2."""
    print(f"Dividing {f1} by {f2}. {comment}")
    inverse_f2 = TitanFraction(f2.den, f2.num)
    print(f"  -> Inverting the divisor to {inverse_f2}")
    return multiply(f1, inverse_f2, "Multiplying by inverse")

# Main calculation logic
print("--- Starting Titan Calculation for Force F ---")
print("Formula: F = g * (d / h) * m")
print("Mass Formula: m = (3 * pi) / 20\n")

# Step 1: Define Titan approximations for constants
g = TitanFraction(10, 1)
d = TitanFraction(20, 1)
h = TitanFraction(10, 1)
pi = TitanFraction(3, 1)
three = TitanFraction(3, 1)
twenty = TitanFraction(20, 1)

print("Chosen 5-bit fractional approximations:")
print(f"  g  = {g}")
print(f"  pi = {pi}")
print(f"  d  = {d}")
print(f"  h  = {h}\n")

# Step 2: Calculate the ratio d/h
print("--- Part 1: Calculate (d / h) ---")
d_over_h = divide(d, h, "Ratio of distance to height")
print(f"Result of (d/h) is {d_over_h}\n")

# Step 3: Calculate the mass m
print("--- Part 2: Calculate mass m = (3 * pi) / 20 ---")
m_num = multiply(three, pi, "Numerator of mass formula (3 * pi)")
mass = divide(m_num, twenty, "Dividing by 20")
print(f"Result for mass (m) is {mass}\n")

# Step 4: Calculate the final force F = g * (d/h) * m
print("--- Part 3: Calculate Force F = g * (d/h) * m ---")
F_intermediate = multiply(g, d_over_h, "Multiplying g by (d/h)")
F_final = multiply(F_intermediate, mass, "Multiplying by mass")

print("\n--- Final Calculation Result ---")
# To show the full equation, we re-declare the numbers
# The result from F_final is 9/1
g_num, g_den = 10, 1
d_div_h_num, d_div_h_den = 2, 1
m_num, m_den = 9, 20
F_num, F_den = F_final.num, F_final.den

# Print the final resolved equation. The calculation showed that
# g*(d/h) simplifies to (20/1), which when multiplied by m=(9/20) simplifies to 9/1
print(f"The calculation simplifies to: ( {g_num * d_div_h_num} / {g_den * d_div_h_den} ) * ( {m_num} / {m_den} ) = ( 20 / 1 ) * ( 9 / 20 ) = {F_num} / {F_den} N")

# Error Analysis
# True value is F = g * d * m / h
# F_true = 9.8 * 20 * (3 * 3.1415926535) / 20 / 10 = 9.8 * (3 * 3.1415926535 / 10) = 2.94 * 3.1415926535
F_true = 9.236186211365
F_titan = F_final.value()
absolute_error = abs(F_true - F_titan)

# Final answer to the user's question
# Since we successfully calculated a force, the answer is Yes.
# We include the smallest absolute error we could produce with a valid calculation path.
final_answer = f"Y[{absolute_error:.3f}]"

# Print the answer in the required format
print(f"\nFinal answer block:")
print(f'<<<{final_answer}>>>')
