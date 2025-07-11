import sympy as sp

def display_magnetic_field_solution():
    """
    This function derives and displays the solution for the magnetic field H.
    The derivation is analytical, and this code serves to format and print the final result.
    """
    # Define symbolic representations of the physical quantities
    # Using strings for clear representation in the final output
    K_0 = "K_0"
    R = "R"
    mu = "μ"
    mu_0 = "μ_0"
    r = "r"
    theta = "θ"
    r_hat = "r̂"
    theta_hat = "θ̂"
    phi_hat = "φ̂"
    z_hat = "ẑ"

    # Based on solving Laplace's equation with magnetostatic boundary conditions,
    # the following expressions for the magnetic field H are obtained.

    # For the region inside the sphere (0 < r < R)
    # The derived field is H_in = (2 * mu_0 * K_0) / (mu + 2 * mu_0) * z_hat
    # To match the format in the answer choices, we rearrange it as:
    # H_in = (2 * mu_0 / mu) * (K_0 / (1 + 2 * mu_0 / mu)) * z_hat
    H_in_coefficient = f"(2 * {mu_0} / {mu}) * ({K_0} / (1 + (2 * {mu_0} / {mu})))"
    H_in_expression = f"H_in(r, θ) = {H_in_coefficient} {z_hat}"

    # For the region outside the sphere (r > R)
    # The derived field is a dipole field.
    # H_out = (K_0 * R^3 / (1 + 2 * mu_0 / mu)) * (1/r^3) * (2*cos(θ)r̂ + sin(θ)θ̂)
    H_out_coefficient = f"({K_0} / (1 + (2 * {mu_0} / {mu})))"
    H_out_spatial_term = f"({R}³ / {r}³)"
    H_out_vector_term = f"(2 * cos({theta}) {r_hat} + sin({theta}) {theta_hat})"
    H_out_expression = f"H_out(r, θ) = {H_out_coefficient} * {H_out_spatial_term} * {H_out_vector_term}"

    # Print the final solution
    print("The magnetic field H(r, θ) inside and outside the sphere is determined as follows:")
    print("\nFor the region inside the sphere (0 < r < R):")
    print(H_in_expression)
    print("\nFor the region outside the sphere (r > R):")
    print(H_out_expression)
    print("\nThese results correspond to Answer Choice E.")


if __name__ == "__main__":
    display_magnetic_field_solution()
