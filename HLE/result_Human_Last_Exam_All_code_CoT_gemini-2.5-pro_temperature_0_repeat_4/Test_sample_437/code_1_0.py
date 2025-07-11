import math

def evaluate_slater_integral(zeta_i, zeta_j):
    """
    Evaluates the integral <phi_i| 1/r |phi_j> for 1s Slater-type orbitals.

    The formula used is: I = (4 * (zeta_i * zeta_j)^(3/2)) / (zeta_i + zeta_j)^2

    Args:
        zeta_i (float): The orbital exponent for the first 1s orbital (phi_i).
        zeta_j (float): The orbital exponent for the second 1s orbital (phi_j).
    """
    print(f"Evaluating the integral for zeta_i = {zeta_i} and zeta_j = {zeta_j}.")

    # Calculate the numerator
    numerator_base = zeta_i * zeta_j
    numerator = 4 * (numerator_base ** 1.5)

    # Calculate the denominator
    denominator_base = zeta_i + zeta_j
    denominator = denominator_base ** 2

    # Calculate the final result
    if denominator == 0:
        # This case is physically not meaningful but good for code robustness
        result = float('inf')
    else:
        result = numerator / denominator

    # Print the calculation step-by-step as requested
    print("\nCalculation:")
    print(f"I = (4 * ({zeta_i} * {zeta_j})^(3/2)) / ({zeta_i} + {zeta_j})^2")
    print(f"I = (4 * {numerator_base}^1.5) / {denominator_base}^2")
    print(f"I = {numerator} / {denominator}")
    print(f"I = {result}")
    
    return result

if __name__ == "__main__":
    # We will evaluate the integral for a 1s orbital of a hydrogen atom.
    # This is the expectation value <1/r>, where phi_i = phi_j.
    # For a hydrogen atom, the orbital exponent zeta is 1.0.
    zeta = 1.0
    
    final_value = evaluate_slater_integral(zeta, zeta)
    
    print(f"\nThe final value of the integral is: {final_value}")
