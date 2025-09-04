import scipy.constants as const
import math

def check_correctness():
    """
    This function checks the correctness of the final answer by recalculating the value
    based on the physical principles described in the problem.
    """
    # --- Setup from the problem ---
    # Physical constants using high-precision values from scipy
    k = const.k  # Coulomb's constant in N m^2 / C^2
    e = const.e  # Elementary charge in C
    
    # Given values from the question
    r = 2.0      # Radius in meters
    q = 2 * e    # Charge of each particle
    
    # For the minimum energy configuration of 12 charges on a sphere (icosahedron),
    # the dimensionless energy constant is E_12. This value is well-established
    # and consistently used in the provided correct analyses.
    E_12 = 49.165
    
    # --- Calculation ---
    # The total potential energy U_total is the sum of two components:
    # 1. Interaction between the central charge and the 12 shell charges: U_cs = 12 * (k * q^2 / r)
    # 2. Mutual interaction between the 12 shell charges: U_ss = E_12 * (k * q^2 / r)
    # The total energy is the sum of these two components.
    # U_total = U_cs + U_ss = (12 + E_12) * (k * q^2 / r)
    
    calculated_energy = (12 + E_12) * (k * q**2 / r)
    
    # --- Verification ---
    # The final answer provided is C, which corresponds to the value 2.822 x 10^-26 J.
    expected_value = 2.822e-26
    
    # The question asks for the answer to be "correct to three decimals".
    # This implies we should check if our calculated value, when formatted similarly, matches.
    # A robust way to compare floating-point numbers is to use a relative tolerance.
    # A tolerance of 0.1% (1e-3) is appropriate for this level of precision, as it
    # accounts for minor rounding differences in constants.
    
    if math.isclose(calculated_energy, expected_value, rel_tol=1e-3):
        return "Correct"
    else:
        # Provide a clear reason for the failure for debugging purposes.
        # Format the calculated value to match the answer's format for easier comparison.
        exponent = math.floor(math.log10(abs(calculated_energy)))
        mantissa = calculated_energy / (10**exponent)
        
        return (f"Incorrect. The calculated energy is approximately {mantissa:.3f} x 10^{exponent} J "
                f"({calculated_energy:.4e} J), which does not match the expected answer of {expected_value:.4e} J "
                f"within the required precision.")

# The final output will be the return value of this function.
print(check_correctness())