import math

def evaluate_1s_integral(zeta_i, zeta_j):
    """
    Evaluates the integral <phi_i | 1/r | phi_j> for 1s Slater orbitals.

    Args:
        zeta_i (float): Orbital exponent for phi_i.
        zeta_j (float): Orbital exponent for phi_j.

    Returns:
        float: The value of the integral.
    """
    numerator = 4 * math.pow(zeta_i * zeta_j, 1.5)
    denominator = math.pow(zeta_i + zeta_j, 2)
    return numerator / denominator

# --- Introduction ---
print("This script evaluates the integral <phi_i | 1/r | phi_j> for 1s Slater-Type Orbitals (STOs).")
print("The general formula for this integral is: I = 4 * (zeta_i * zeta_j)^(3/2) / (zeta_i + zeta_j)^2\n")

# --- General Case Calculation ---
print("--- Case 1: General case with different orbitals ---")
zeta_i_val = 1.0
zeta_j_val = 1.2
print(f"Let's use example values: zeta_i = {zeta_i_val}, zeta_j = {zeta_j_val}\n")

# Step-by-step calculation printout
print("Calculation steps:")
print(f"I = 4 * ({zeta_i_val} * {zeta_j_val})^(3/2) / ({zeta_i_val} + {zeta_j_val})^2")

step1_num = zeta_i_val * zeta_j_val
step1_den = zeta_i_val + zeta_j_val
print(f"I = 4 * ({step1_num})^(1.5) / ({step1_den})^2")

step2_num = math.pow(step1_num, 1.5)
step2_den = math.pow(step1_den, 2)
print(f"I = 4 * {step2_num:.4f} / {step2_den:.4f}")

step3_num = 4 * step2_num
print(f"I = {step3_num:.4f} / {step2_den:.4f}")

final_result_1 = evaluate_1s_integral(zeta_i_val, zeta_j_val)
print(f"Final Result: I = {final_result_1:.6f}")
print("-" * 50)

# --- Special Case Calculation ---
print("\n--- Case 2: Special case for the expectation value <phi | 1/r | phi> ---")
print("In this case, phi_i = phi_j, which means zeta_i = zeta_j = zeta.")
print("The formula simplifies beautifully to just: I = zeta.\n")

zeta_val = 1.5
print(f"Let's use an example value: zeta = {zeta_val}")
print(f"The expectation value is simply the value of zeta.")
print(f"Final Result: <1/r> = {zeta_val}\n")

print("Let's verify this special case using the general formula function:")
verification_result = evaluate_1s_integral(zeta_val, zeta_val)
print(f"Calculation with general formula for zeta = {zeta_val}:")
print(f"I = 4 * ({zeta_val} * {zeta_val})^(1.5) / ({zeta_val} + {zeta_val})^2 = {verification_result}")
print("The result matches, confirming the simplification is correct.")
