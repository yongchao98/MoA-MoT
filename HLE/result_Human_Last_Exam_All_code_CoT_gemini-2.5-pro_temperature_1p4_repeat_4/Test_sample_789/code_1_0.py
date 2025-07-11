import math

# Helper function for greatest common divisor, needed for simplification
def gcd(a, b):
    while b:
        a, b = b, a % b
    return a

# Represents a fraction with 5-bit integer components
class TitanFraction:
    def __init__(self, num, den):
        if not (0 <= num <= 31 and 1 <= den <= 31):
            # This check is for initial setup. Overflows are handled in multiply.
            raise ValueError(f"Numerator {num} or denominator {den} out of 5-bit range (0-31)")
        self.num = num
        self.den = den

    def __str__(self):
        return f"{self.num}/{self.den}"

    def value(self):
        return self.num / self.den

# Simulates the Titan multiplication with cross-simplification
def multiply(f1, f2):
    print(f"Multiplying {f1} * {f2}")
    
    # Cross-simplify f1.num and f2.den
    common_divisor_1 = gcd(f1.num, f2.den)
    num1 = f1.num // common_divisor_1
    den2 = f2.den // common_divisor_1
    print(f"Cross-simplifying {f1.num} and {f2.den} by {common_divisor_1} -> {num1} and {den2}")

    # Cross-simplify f2.num and f1.den
    common_divisor_2 = gcd(f2.num, f1.den)
    num2 = f2.num // common_divisor_2
    den1 = f1.den // common_divisor_2
    print(f"Cross-simplifying {f2.num} and {f1.den} by {common_divisor_2} -> {num2} and {den1}")

    # Multiply the simplified components
    res_num = num1 * num2
    res_den = den1 * den2
    print(f"Resulting multiplication: ({num1}*{num2}) / ({den1}*{den2}) = {res_num}/{res_den}")

    # Check for overflow in the final result before simplification
    if res_num > 31 or res_den > 31:
        print(f"Error: Overflow detected! Result {res_num}/{res_den} is not representable.")
        return None
    
    # Final simplification of the result itself (e.g., 12/10 -> 6/5)
    final_common_divisor = gcd(res_num, res_den)
    if final_common_divisor > 1:
      print(f"Simplifying final fraction {res_num}/{res_den} by {final_common_divisor}")
      res_num //= final_common_divisor
      res_den //= final_common_divisor

    print(f"Final result of multiplication: {res_num}/{res_den}\n")
    return TitanFraction(res_num, res_den)

def calculate_mass():
    print("Calculating mass of the rock.")
    print("Formula: mass = ρ * (4/3) * π * r³")
    
    # Step 1: Define the components as 5-bit fractions
    # ρ = 0.9 = 9/10
    # r = 0.5 = 1/2 -> r³ = 1/8
    # π is approximated as 28/9 to minimize error while avoiding overflow
    rho = TitanFraction(9, 10)
    four_thirds = TitanFraction(4, 3)
    pi_approx = TitanFraction(28, 9)
    r_cubed = TitanFraction(1, 8)
    
    print(f"Using ρ = {rho}, r³ = {r_cubed}, and π ≈ {pi_approx}\n")
    
    # Step 2: Perform calculations in an order that avoids overflows
    # Path: mass = ( (9/10) * (4/3) ) * ( (28/9) * (1/8) )
    
    print("--- Calculating first part: (ρ * 4/3) ---")
    term1 = multiply(rho, four_thirds)
    
    print("--- Calculating second part: (π * r³) ---")
    term2 = multiply(pi_approx, r_cubed)
    
    print("--- Calculating final mass: (part1 * part2) ---")
    final_mass = multiply(term1, term2)
    
    print("--------------------------------------------------")
    print("Final Calculation Path:")
    print(f"mass = ( ({rho.num}/{rho.den}) * ({four_thirds.num}/{four_thirds.den}) ) * ( ({pi_approx.num}/{pi_approx.den}) * ({r_cubed.num}/{r_cubed.den}) )")
    print(f"mass = ( {term1.num}/{term1.den} ) * ( {term2.num}/{term2.den} )")
    print(f"Final mass = {final_mass.num} / {final_mass.den}")
    
    # Step 3: Calculate the absolute error
    true_mass_val = 0.9 * (4/3) * math.pi * (0.5**3)
    calculated_mass_val = final_mass.value()
    abs_error = abs(true_mass_val - calculated_mass_val)
    
    # The question asks for the value of e, rounded to 0.001
    rounded_error = round(abs_error, 3)
    
    # The final output needs to be in a special format <<<content>>>
    # This print is for the user to see the value.
    print(f"\nSmallest absolute error (e) = {abs_error:.7f}")
    print(f"Rounded to 0.001, e = {rounded_error}")
    
    return rounded_error

# Run the calculation
final_answer = calculate_mass()

# The final return as requested by the format
print(f"\n<<<e={final_answer}>>>")
