import math

def solve_magnetic_shell():
    """
    Calculates the required permeability of a cylindrical magnetic shell to cloak
    an interior diamagnetic core, and the resulting magnetic field inside the core.
    """

    # --- User-defined parameters for the calculation ---
    # Inner radius of the shell in meters
    R1 = 0.1
    # Outer radius of the shell in meters
    R2 = 0.2
    # Magnitude of the uniform external magnetic field in A/m
    H0 = 1.0

    # --- Physical Constants ---
    # Permeability of free space in H/m
    mu_0 = 4 * math.pi * 1e-7

    # --- Introduction and Assumptions ---
    print("This solution determines the properties of a cylindrical shell that 'cloaks' an interior object from an external magnetic field.")
    print("The condition that the external field is not distorted, while requiring a non-trivial solution (mu != mu_0),")
    print("implies that the interior region (rho < R1) is not a vacuum. We assume it is a perfect diamagnet (e.g., a superconductor) with permeability mu_1 = 0.\n")
    print(f"Given parameters: R1 = {R1} m, R2 = {R2} m, H0 = {H0} A/m.")

    # --- Part 1: Calculate the required permeability of the shell (mu) ---
    print("\n" + "="*60)
    print("1. Required Permeability of the Shell (mu)")
    print("="*60)
    print("The required permeability 'mu' is given by the formula:")
    print("mu = mu_0 * (R2^2 + R1^2) / (R2^2 - R1^2)\n")

    # Calculation with numbers
    r1_sq = R1**2
    r2_sq = R2**2
    numerator_mu = r2_sq + r1_sq
    denominator_mu = r2_sq - r1_sq
    
    if denominator_mu <= 0:
        print("Error: R2 must be greater than R1.")
        return

    mu = mu_0 * numerator_mu / denominator_mu

    print("Plugging in the values:")
    print(f"mu = ({mu_0:.5e}) * (({R2})^2 + ({R1})^2) / (({R2})^2 - ({R1})^2)")
    print(f"mu = ({mu_0:.5e}) * ({r2_sq:.4f} + {r1_sq:.4f}) / ({r2_sq:.4f} - {r1_sq:.4f})")
    print(f"mu = ({mu_0:.5e}) * ({numerator_mu:.4f}) / ({denominator_mu:.4f})")
    print(f"\nResult: mu = {mu:.5e} H/m")

    # --- Part 2: Calculate the magnetic field in the interior (H_int) ---
    print("\n" + "="*60)
    print("2. Magnetic Field in the Interior (H_int)")
    print("="*60)
    print("The magnetic field in the interior is uniform, points in the x-direction,")
    print("and its magnitude |H_int| is given by:")
    print("|H_int| = H0 * (2 * R2^2) / (R1^2 + R2^2)\n")

    # Calculation with numbers
    numerator_h = 2 * r2_sq
    denominator_h = r1_sq + r2_sq
    h_int_mag = H0 * numerator_h / denominator_h

    print("Plugging in the values:")
    print(f"|H_int| = {H0} * (2 * ({R2})^2) / (({R1})^2 + ({R2})^2)")
    print(f"|H_int| = {H0} * (2 * {r2_sq:.4f}) / ({r1_sq:.4f} + {r2_sq:.4f})")
    print(f"|H_int| = {H0} * ({numerator_h:.4f}) / ({denominator_h:.4f})")
    print(f"\nResult: |H_int| = {h_int_mag:.4f} A/m")

    print("\n" + "="*60)
    print("Final Answer Summary:")
    print(f"The required permeability is mu = {mu:.5e} H/m.")
    print(f"The interior magnetic field is H_int = {h_int_mag:.4f} x_hat A/m.")
    print("="*60)

# Execute the function to solve the problem
solve_magnetic_shell()