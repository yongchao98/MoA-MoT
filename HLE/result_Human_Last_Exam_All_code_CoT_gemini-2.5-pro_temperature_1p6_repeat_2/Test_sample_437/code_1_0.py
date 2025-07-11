import math

def evaluate_integral(zeta_i, zeta_j):
    """
    Evaluates the integral <phi_i | 1/r | phi_j> for 1s Slater orbitals.
    The analytical formula is: 4 * (zeta_i * zeta_j)^(3/2) / (zeta_i + zeta_j)^2
    """
    
    # Calculate numerator
    zeta_product = zeta_i * zeta_j
    numerator = 4 * (zeta_product ** 1.5)
    
    # Calculate denominator
    zeta_sum = zeta_i + zeta_j
    denominator = zeta_sum ** 2
    
    # Calculate final result
    result = numerator / denominator
    
    # Print the detailed calculation steps
    print(f"Evaluating the integral <phi|1/r|phi> for a 1s Slater orbital.")
    print(f"We are considering the case where the orbital exponent zeta = {zeta_i}\n")
    print("The general formula for the integral <phi_i|1/r|phi_j> is:")
    print("I = 4 * (zeta_i * zeta_j)^(3/2) / (zeta_i + zeta_j)^2\n")

    print("--- Calculation Steps ---")
    print(f"1. Set exponents for the expectation value calculation:")
    print(f"   zeta_i = {zeta_i}")
    print(f"   zeta_j = {zeta_j}\n")
    
    print(f"2. Calculate the numerator:")
    print(f"   Numerator = 4 * ({zeta_i} * {zeta_j})^(3/2)")
    print(f"             = 4 * ({zeta_product})^(1.5)")
    print(f"             = {numerator}\n")
    
    print(f"3. Calculate the denominator:")
    print(f"   Denominator = ({zeta_i} + {zeta_j})^2")
    print(f"               = ({zeta_sum})^2")
    print(f"               = {denominator}\n")

    print(f"4. Calculate the final result:")
    print(f"   Result = Numerator / Denominator")
    print(f"          = {numerator} / {denominator}")
    print(f"          = {result}\n")
    
    print(f"As derived analytically, the result {result} is equal to the orbital exponent zeta = {zeta_i}.")

    return result

# We evaluate the expectation value for a 1s orbital.
# For the hydrogen atom, the 1s orbital has an exponent zeta = 1.0.
zeta = 1.0
final_value = evaluate_integral(zeta, zeta)