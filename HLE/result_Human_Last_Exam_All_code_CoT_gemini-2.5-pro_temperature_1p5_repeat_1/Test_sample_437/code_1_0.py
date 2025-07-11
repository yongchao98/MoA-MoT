import math

def evaluate_integral(zeta_i, zeta_j):
    """
    Evaluates the integral <phi_i| 1/r |phi_j> for 1s Slater orbitals.
    The general formula is: 4 * (zeta_i * zeta_j)^(3/2) / (zeta_i + zeta_j)^2
    """
    
    print("--- Calculation Steps ---")

    # Step 1: Calculate intermediate terms
    zi_x_zj = zeta_i * zeta_j
    zi_p_zj = zeta_i + zeta_j
    print(f"1. Let zeta_i = {zeta_i} and zeta_j = {zeta_j}")
    print(f"2. The product of exponents: zeta_i * zeta_j = {zeta_i} * {zeta_j} = {zi_x_zj}")
    print(f"3. The sum of exponents: zeta_i + zeta_j = {zeta_i} + {zeta_j} = {zi_p_zj}")

    # Step 2: Calculate the numerator
    num_base = zi_x_zj
    num_exponent = 1.5  # 3/2
    num_base_pow = math.pow(num_base, num_exponent)
    numerator = 4 * num_base_pow
    print(f"4. Numerator calculation: 4 * ({num_base})^(3/2) = 4 * {num_base_pow:.4f} = {numerator:.4f}")

    # Step 3: Calculate the denominator
    denominator = math.pow(zi_p_zj, 2)
    print(f"5. Denominator calculation: ({zi_p_zj})^2 = {denominator:.4f}")

    # Step 4: Final division
    result = numerator / denominator
    print(f"6. Final result: {numerator:.4f} / {denominator:.4f} = {result:.4f}")
    
    return result

# The prompt is for "a 1s Slater orbital", meaning phi_i = phi_j.
# This implies their exponents are equal, zeta_i = zeta_j = zeta.
# For a Hydrogen atom's 1s orbital, a simple choice is zeta = 1.0.
zeta = 1.0

print(f"Evaluating the integral for a single 1s Slater orbital with exponent zeta = {zeta}")
print("This corresponds to the case where phi_i = phi_j, so zeta_i = zeta_j = zeta.\n")

final_value = evaluate_integral(zeta, zeta)

print("\n--- Conclusion ---")
print(f"The calculated value of the integral is {final_value:.4f}.")
print(f"This matches the analytical result, which is the value of the exponent zeta ({zeta}).")
