import math

def evaluate_1s_slater_integral(zeta_i, zeta_j):
    """
    Evaluates the integral <phi_i| 1/r |phi_j> for 1s Slater-Type Orbitals (STOs).

    Args:
        zeta_i (float): The orbital exponent for the first 1s STO (phi_i).
        zeta_j (float): The orbital exponent for the second 1s STO (phi_j).

    Returns:
        float: The calculated value of the integral.
    """
    print("Evaluating the integral <phi_i| 1/r |phi_j> for 1s Slater-Type Orbitals.")
    print(f"The orbital exponents are zeta_i = {zeta_i} and zeta_j = {zeta_j}.")
    print("-" * 50)

    # The analytical formula for the integral is:
    # I = 4 * (zeta_i * zeta_j)^(3/2) / (zeta_i + zeta_j)^2
    
    # Calculate intermediate parts of the formula
    numerator_base = zeta_i * zeta_j
    denominator_base = zeta_i + zeta_j
    
    numerator = 4 * (numerator_base ** 1.5)
    denominator = denominator_base ** 2
    
    result = numerator / denominator

    # Print the equation with the numbers plugged in
    print("The formula for the integral is:")
    print("I = 4 * (zeta_i * zeta_j)^(3/2) / (zeta_i + zeta_j)^2\n")
    
    print("Substituting the given values:")
    print(f"I = 4 * ({zeta_i} * {zeta_j})^(3/2) / ({zeta_i} + {zeta_j})^2")
    print(f"I = 4 * ({numerator_base})^(1.5) / ({denominator_base})^2")
    print(f"I = {numerator} / {denominator}")
    print("-" * 50)
    print(f"The final result of the integral is: {result}")
    
    return result

# --- Example Calculation ---
# We will evaluate the expectation value <1s| 1/r |1s> for a Helium atom,
# using an effective orbital exponent from variational theory, zeta = 1.6875.
# In this case, zeta_i = zeta_j.
zeta_val = 1.6875

# Call the function with the chosen zeta values
final_answer = evaluate_1s_slater_integral(zeta_val, zeta_val)

# The final answer is also printed inside the function.
# For the special case where zeta_i = zeta_j = zeta, the result is simply zeta.
# Let's verify: 1.6875
# The code should output this value.