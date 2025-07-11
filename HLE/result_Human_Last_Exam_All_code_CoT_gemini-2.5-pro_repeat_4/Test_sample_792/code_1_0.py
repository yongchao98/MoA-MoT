import math

def gcd(a, b):
    """Computes the greatest common divisor of a and b."""
    return math.gcd(a, b)

def simplify_fraction(n, d):
    """Simplifies a fraction by dividing numerator and denominator by their GCD."""
    if d == 0:
        raise ValueError("Denominator cannot be zero.")
    common_divisor = gcd(n, d)
    return n // common_divisor, d // common_divisor

def check_5bit(num):
    """Checks if a number is a valid 5-bit integer (0-31)."""
    if not (0 <= num <= 31):
        raise ValueError(f"Constraint Violated: {num} is not a valid 5-bit integer.")

def multiply_fractions(f1_n, f1_d, f2_n, f2_d):
    """
    Multiplies two fractions, simplifying before computing the final product
    to avoid violating 5-bit constraints on intermediate products.
    """
    # Simplify by cross-cancellation first
    g1 = gcd(f1_n, f2_d)
    f1_n //= g1
    f2_d //= g1
    
    g2 = gcd(f2_n, f1_d)
    f2_n //= g2
    f1_d //= g2

    # Now compute the product of the simplified terms
    res_n = f1_n * f2_n
    res_d = f1_d * f2_d

    # Check if the immediate result violates constraints
    check_5bit(res_n)
    check_5bit(res_d)

    # The result is already simplified as much as possible by cross-cancellation
    return res_n, res_d


# --- Titan Calculation ---

print("Starting Titan computation for the force F.")
print("The derived formula is: F = (3/10) * g * pi")
print("-" * 30)

# Step 1: Define initial fractional constants
# Base term from the formula F = (3/10) * g * pi
f_base_n, f_base_d = 3, 10
# Approximation for g (9.8 m/s^2) -> 28/3 = 9.33...
g_n, g_d = 28, 3
# Approximation for pi (~3.14) -> 10/3 = 3.33...
pi_n, pi_d = 10, 3

print(f"Using g = {g_n}/{g_d} and pi = {pi_n}/{pi_d}")

# Step 2: Calculate the first part of the equation: (3/10) * (28/3)
print(f"\nStep 1: Calculate ({f_base_n}/{f_base_d}) * ({g_n}/{g_d})")
try:
    # Manual simplification for clarity
    temp_n_display = f_base_n * g_n  # 84
    temp_d_display = f_base_d * g_d  # 30
    
    # Simplify by cancelling common factor 3
    simplified_n1, simplified_d1 = 28, 10
    print(f"  Intermediate product: ({temp_n_display}/{temp_d_display}). This simplifies by cancelling '3' to ({simplified_n1}/{simplified_d1}).")

    # Further simplify 28/10 to 14/5
    interim_n, interim_d = simplify_fraction(simplified_n1, simplified_d1)
    check_5bit(interim_n)
    check_5bit(interim_d)
    print(f"  Further simplification yields: {interim_n}/{interim_d}. This is a valid intermediate result.")

except ValueError as e:
    print(f"  Calculation failed: {e}")
    exit()

# Step 3: Multiply the intermediate result by the fraction for pi
print(f"\nStep 2: Calculate ({interim_n}/{interim_d}) * ({pi_n}/{pi_d})")
try:
    # Manual simplification for clarity
    temp_n_display = interim_n * pi_n # 14 * 10 = 140
    temp_d_display = interim_d * pi_d # 5 * 3 = 15

    # Simplify by cancelling common factor 5
    final_n = 14 * 2 # 28
    final_d = 3      # 3
    print(f"  Intermediate product: ({temp_n_display}/{temp_d_display}). This simplifies by cancelling '5' to ({final_n}/{final_d}).")

    check_5bit(final_n)
    check_5bit(final_d)
    print(f"  Final result is {final_n}/{final_d}. This is a valid Titan fraction.")

except ValueError as e:
    print(f"  Calculation failed: {e}")
    exit()

print("-" * 30)
print("Final Calculation Summary:")
# The problem asks to output each number in the final equation.
print(f"F = ({f_base_n}/{f_base_d}) * ({g_n}/{g_d}) * ({pi_n}/{pi_d}) = {final_n}/{final_d}")

# Final error calculation
f_true = 2.94 * math.pi
f_titan = final_n / final_d
error = abs(f_titan - f_true)
print(f"\nTrue Force Value: ~{f_true:.4f} N")
print(f"Titan Calculated Force: {final_n}/{final_d} = {f_titan:.4f} N")
print(f"Absolute Error: {error:.4f}")

# Check if the rock hits the lion
# Target range for force F is [8.774, 9.698] N
f_min_true = (19/20) * f_true
f_max_true = (21/20) * f_true
if f_min_true <= f_titan <= f_max_true:
    print("\nThe calculated force is within the required range to hit the lion.")
else:
    print("\nThe calculated force is NOT within the required range to hit the lion.")
