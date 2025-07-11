def solve_magnetic_shell():
    """
    This function calculates and prints the symbolic formulas for the required
    permeability of a cylindrical shell and the resulting interior magnetic field
    such that an external uniform magnetic field is not distorted.
    """

    # Define symbolic variables as strings for printing formulas
    R1_str = "R1"
    R2_str = "R2"
    H0_str = "H0"
    mu0_str = "mu_0"

    # The derivation leads to two possible solutions for the relative permeability (mu/mu_0).
    # We present the solution where mu > mu_0, which is common in such problems.
    
    # --- Permeability Calculation ---
    # The required relative permeability mu_r = mu / mu_0 is given by:
    # mu_r = (R2 + R1) / (R2 - R1)
    # Therefore, the absolute permeability mu is:
    # mu = mu_0 * (R2 + R1) / (R2 - R1)
    
    print("The required permeability of the shell material is given by the formula:")
    # The instruction "output each number in the final equation" is interpreted
    # as showing all coefficients, which are 1 in this case.
    print(f"mu = {mu0_str} * ( (1 * {R2_str} + 1 * {R1_str}) / (1 * {R2_str} - 1 * {R1_str}) )")
    print("-" * 20)

    # --- Interior Magnetic Field Calculation ---
    # For this permeability, the magnetic field inside the shell (rho < R1)
    # is found to be uniform and is given by:
    # H_int = H0 * (R2 / R1) * x_hat
    
    print("The corresponding magnetic field in the interior region (rho < R1) is:")
    print(f"H_int = {H0_str} * ({R2_str} / {R1_str}) * x_hat")
    print("\nNote: 'x_hat' is the unit vector in the x-direction.")

solve_magnetic_shell()

# The final answer consists of the two derived formulas.
# Permeability: mu = mu_0 * (R2 + R1) / (R2 - R1)
# Interior Field: H_int = H0 * (R2 / R1) * x_hat
# Since the output is symbolic, we format it as a string.
final_answer_str = "mu = mu_0 * (R2 + R1) / (R2 - R1) and H_int = H0 * (R2 / R1) * x_hat"
print(f"\n<<<Final Answer Summary: {final_answer_str}>>>")