import numpy as np
from scipy.integrate import quad

def evaluate_1s_nuclear_attraction(zeta_i, zeta_j):
    """
    Evaluates the nuclear attraction integral <phi_i| 1/r |phi_j> for two
    1s Slater-Type Orbitals (STOs) centered on the same nucleus.

    Args:
        zeta_i (float): The orbital exponent of the first 1s orbital (phi_i).
        zeta_j (float): The orbital exponent of the second 1s orbital (phi_j).

    The integral is defined as:
    I = Integral( phi_i*(r) * (1/r) * phi_j(r) * d_tau ) over all space.

    The normalized 1s STO is: phi(r) = (zeta^3 / pi)^(1/2) * exp(-zeta*r)

    The integral can be solved analytically, yielding the formula:
    I = 4 * (zeta_i * zeta_j)^(3/2) / (zeta_i + zeta_j)^2

    This function calculates the result using the analytical formula and also
    verifies it by performing numerical integration.
    """

    # --- Analytical Calculation ---
    numerator = 4 * (zeta_i * zeta_j)**1.5
    denominator = (zeta_i + zeta_j)**2
    analytical_result = numerator / denominator

    # --- Numerical Calculation for Verification ---
    # The full integrand in spherical coordinates, after integrating over angles, is:
    # Integrand(r) = phi_i(r) * (1/r) * phi_j(r) * 4*pi*r^2
    # phi_i(r) = (zeta_i^3 / pi)^0.5 * exp(-zeta_i*r)
    # phi_j(r) = (zeta_j^3 / pi)^0.5 * exp(-zeta_j*r)
    # Substituting these in simplifies the integrand to:
    # Integrand(r) = 4 * (zeta_i*zeta_j)^1.5 * r * exp(-(zeta_i+zeta_j)*r)
    
    radial_integrand = lambda r: 4 * (zeta_i * zeta_j)**1.5 * r * np.exp(-(zeta_i + zeta_j) * r)
    numerical_result, error = quad(radial_integrand, 0, np.inf)

    # --- Output the results ---
    print(f"Evaluating the integral <phi_i| 1/r |phi_j> for 1s STOs with exponents:")
    print(f"zeta_i = {zeta_i}")
    print(f"zeta_j = {zeta_j}")
    print("-" * 50)
    
    print("ANALYTICAL RESULT:")
    print("The final equation for the integral is:")
    # Printing each number in the final equation as requested
    print(f"I = {4.0} * ({zeta_i} * {zeta_j})**{1.5} / ({zeta_i} + {zeta_j})**{2.0}")
    print(f"Analytical Value = {analytical_result}")
    print("-" * 50)

    print("NUMERICAL VERIFICATION:")
    print("Numerical Value = {}".format(numerical_result))
    print("Absolute Error in Numerical Integration = {}".format(error))
    print("-" * 50)
    
    # In the special case where zeta_i = zeta_j = zeta, the formula simplifies to zeta.
    if zeta_i == zeta_j:
        print(f"Note: For the case zeta_i = zeta_j = {zeta_i}, the integral simplifies to just zeta, which is {analytical_result}.")

# --- Example Usage ---
# You can change these values to evaluate the integral for different orbitals.
# Case 1: Different exponents
print("### General Case (zeta_i != zeta_j) ###")
zeta_i_val_1 = 2.0
zeta_j_val_1 = 1.5
result_1 = evaluate_1s_nuclear_attraction(zeta_i_val_1, zeta_j_val_1)
print("\n")

# Case 2: Identical exponents (e.g., for the expectation value <phi|1/r|phi>)
print("### Special Case (zeta_i = zeta_j) ###")
# For a Hydrogen atom, zeta is 1.0. For a He atom, a common value is 1.6875.
zeta_val_2 = 1.6875
result_2 = evaluate_1s_nuclear_attraction(zeta_val_2, zeta_val_2)
