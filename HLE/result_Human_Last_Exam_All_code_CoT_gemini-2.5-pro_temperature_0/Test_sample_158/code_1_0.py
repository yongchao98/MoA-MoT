import math

def solve_magnetic_shell():
    """
    Calculates the required permeability and interior magnetic field for a
    cylindrical shell that does not distort an external magnetic field.
    """
    # Define the parameters of the cylindrical shell.
    # Using example values for demonstration purposes.
    R1 = 1.0  # Inner radius
    R2 = 2.0  # Outer radius
    H0 = 1.0  # Magnitude of the applied uniform magnetic field

    print(f"Solving for a cylindrical shell with R1 = {R1} and R2 = {R2}.")
    print(f"The applied external field is H0 = {H0} in the x-direction.\n")

    # The condition that the external field is not distorted leads to the
    # following equation for the relative permeability mu_r = mu / mu_0:
    # (mu_r - 1)^2 = (R2/R1)^2 * (mu_r + 1)^2
    # This equation has two non-trivial solutions.

    k = R2 / R1

    # --- Solution 1 ---
    # This solution gives a negative relative permeability with magnitude < 1.
    mu_r1 = -(k - 1) / (k + 1)
    # The corresponding interior field is H_int = -H0 * (R2/R1) * x_hat
    H_int_factor1 = -k

    print("--- Solution 1 ---")
    print("The required relative permeability mu_r = mu/mu_0 is:")
    print(f"mu_r = ({R1} - {R2}) / ({R1} + {R2}) = {mu_r1:.4f}")
    print("\nThe resulting magnetic field in the interior region (rho < R1) is:")
    print(f"H_int = - H0 * (R2/R1) * x_hat")
    print(f"H_int = -{H0:.1f} * ({R2:.1f}/{R1:.1f}) * x_hat = {H_int_factor1:.4f} * H0 * x_hat\n")

    # --- Solution 2 ---
    # This solution gives a negative relative permeability with magnitude > 1.
    mu_r2 = -(k + 1) / (k - 1)
    # The corresponding interior field is H_int = H0 * (R2/R1) * x_hat
    H_int_factor2 = k

    print("--- Solution 2 ---")
    print("The required relative permeability mu_r = mu/mu_0 is:")
    print(f"mu_r = -({R1} + {R2}) / ({R2} - {R1}) = {mu_r2:.4f}")
    print("\nThe resulting magnetic field in the interior region (rho < R1) is:")
    print(f"H_int = H0 * (R2/R1) * x_hat")
    print(f"H_int = {H0:.1f} * ({R2:.1f}/{R1:.1f}) * x_hat = {H_int_factor2:.4f} * H0 * x_hat")

# Execute the function
solve_magnetic_shell()