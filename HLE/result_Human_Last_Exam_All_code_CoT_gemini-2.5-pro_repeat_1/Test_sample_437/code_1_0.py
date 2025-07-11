import math

def evaluate_1s_integral():
    """
    Evaluates the integral <phi_i | 1/r | phi_j> for 1s Slater orbitals
    on the same center for a specific case and prints the steps.
    """
    # Define the orbital exponents for the example case.
    # These could represent two different 1s basis functions in a calculation.
    zeta_i = 1.0
    zeta_j = 2.0

    print(f"This script evaluates the integral I = <phi_i| 1/r |phi_j> for 1s Slater orbitals.")
    print(f"We use the analytical formula: I = 4 * (zeta_i * zeta_j)**(3/2) / (zeta_i + zeta_j)**2")
    print("\n--- Calculation for a Specific Case ---")
    print(f"Let's evaluate the integral for zeta_i = {zeta_i} and zeta_j = {zeta_j}.")

    # Calculate the components of the formula
    zeta_product = zeta_i * zeta_j
    zeta_sum = zeta_i + zeta_j

    numerator = 4 * (zeta_product)**1.5
    denominator = zeta_sum**2
    
    integral_value = numerator / denominator

    # Print the final equation with all intermediate numbers
    print("\n--- Final Equation ---")
    print(f"I = (4 * ({zeta_i} * {zeta_j})**1.5) / ({zeta_i} + {zeta_j})**2")
    print(f"I = (4 * {zeta_product}**1.5) / {zeta_sum}**2")
    print(f"I = {numerator} / {denominator}")
    print(f"I = {integral_value}")

if __name__ == "__main__":
    evaluate_1s_integral()