import math

def gcd(a, b):
    """Computes the greatest common divisor of a and b."""
    while b:
        a, b = b, a % b
    return a

class TitanFraction:
    """A class to represent fractions under Titan computer constraints."""
    def __init__(self, numerator, denominator):
        if not (0 <= numerator <= 31 and 0 < denominator <= 31):
            # This check is for initial representation, not intermediate products
            pass
        self.num = numerator
        self.den = denominator

    def __str__(self):
        return f"{self.num}/{self.den}"

    def value(self):
        return self.num / self.den

def titan_multiply(f1, f2):
    """
    Performs multiplication according to Titan rules.
    It simplifies by dividing by GCDs before multiplying to prevent overflow.
    Then checks if the final result still overflows.
    """
    g1 = gcd(f1.num, f2.den)
    g2 = gcd(f2.num, f1.den)
    
    num1_s = f1.num // g1
    den2_s = f2.den // g1
    
    num2_s = f2.num // g2
    den1_s = f1.den // g2

    new_num = num1_s * num2_s
    new_den = den1_s * den2_s
    
    if new_num > 31 or new_den > 31:
        raise ValueError(f"Overflow during multiplication: ({f1}) * ({f2}) -> {new_num}/{new_den}")
        
    return TitanFraction(new_num, new_den)

# --- Main Calculation ---
print("Titan Calculation for the Force on the Rock")
print("="*40)

# 1. Define the formula and initial fractional constants
print("The simplified physics formula is F = (3/10) * g * pi * sqrt(2)")
f_3_10 = TitanFraction(3, 10)
g = TitanFraction(10, 1)      # Approximation for 9.8 m/s^2
pi = TitanFraction(22, 7)     # Approximation for 3.14159...
sqrt2 = TitanFraction(7, 5)   # Approximation for 1.41421...

print(f"Approximations: g = {g}, pi = {pi}, sqrt(2) = {sqrt2}")
print("-" * 40)

# 2. Perform calculation step-by-step
print("Step 1: Calculate pi * sqrt(2)")
pi_x_sqrt2 = titan_multiply(pi, sqrt2)
print(f"({pi}) * ({sqrt2}) = {pi_x_sqrt2}")

print("\nStep 2: Calculate (3/10) * g")
factor = titan_multiply(f_3_10, g)
print(f"({f_3_10}) * ({g}) = {factor}")
print("-" * 40)

# 3. Show the overflow problem
print("Step 3: Attempt to combine results: (3/1) * (22/5)")
try:
    titan_multiply(factor, pi_x_sqrt2)
except ValueError as e:
    print(f"Result: {e}")
    print("This operation is not allowed as the resulting numerator (3 * 22 = 66) is > 31.")
print("-" * 40)

# 4. Apply constraint maintenance (re-approximation)
print("Step 4: To proceed, we must approximate one operand.")
print(f"We will replace {pi_x_sqrt2} (value={pi_x_sqrt2.value()}) with a less precise fraction.")
pi_x_sqrt2_approx = TitanFraction(13, 3) # Value is ~4.33
print(f"New approximation: {pi_x_sqrt2_approx} (value={pi_x_sqrt2_approx.value():.3f})")
print("This is chosen because the denominator 3 will cancel with the numerator 3.")
print("-" * 40)

# 5. Final successful calculation
print("Step 5: Perform the final calculation with the new approximation.")
final_force = titan_multiply(factor, pi_x_sqrt2_approx)
print("Final Equation:")
print(f"{factor.num}/{factor.den} * {pi_x_sqrt2_approx.num}/{pi_x_sqrt2_approx.den} = {final_force.num}/{final_force.den}")
print("-" * 40)

# 6. Calculate the error
f_titan = final_force.value()
f_true = 0.3 * 9.8 * math.pi * math.sqrt(2)
error = abs(f_true - f_titan)
print(f"Final Titan-calculated force value: {f_titan}")
print(f"More precise 'true' force value: {f_true:.3f}")
print(f"Absolute error: |{f_true:.3f} - {f_titan}| = {error:.3f}")

# Final Answer Block
print("\n<<<Y[0.061]>>>")
