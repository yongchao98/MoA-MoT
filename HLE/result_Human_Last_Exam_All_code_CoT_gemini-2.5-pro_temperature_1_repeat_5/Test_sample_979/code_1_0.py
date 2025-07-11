def solve_magnetic_field():
    """
    This function presents the derived solution for the magnetic field H
    inside and outside a spherical shell with a given surface current.
    The derivation involves solving Laplace's equation for the magnetic scalar potential
    and applying the boundary conditions at the sphere's surface.
    """
    # Define symbolic representations for the physical quantities
    K0 = "K_0"
    mu = "μ"
    mu0 = "μ_0"
    R = "R"
    r = "r"
    theta = "θ"
    z_hat = "ẑ"
    r_hat = "r̂"
    theta_hat = "θ̂"

    # The problem is solved by finding coefficients for the magnetic scalar potential.
    # The boundary conditions at r=R are:
    # 1) μ * H_in,r = μ_0 * H_out,r
    # 2) H_out,θ - H_in,θ = K_0 * sin(θ)
    # Solving this system for the l=1 mode gives the coefficients and the fields.
    
    # Expression for the H-field inside the sphere (r < R)
    # H_in = ( (2 * μ_0 * K_0) / (μ + 2 * μ_0) ) * ẑ
    # Rearranging to match the answer choices format:
    # H_in = (2 * μ_0 / μ) * K_0 / (1 + 2 * μ_0 / μ) * ẑ
    H_in_coeff = f"(2*{mu0}/{mu}) * {K0} / (1 + (2*{mu0}/{mu}))"
    H_in_vector = f"{H_in_coeff} {z_hat}"

    # Expression for the H-field outside the sphere (r > R)
    # H_out = (μ * K_0 * R^3 / (μ + 2 * μ_0)) * (1/r^3) * (2*cos(θ)r̂ + sin(θ)θ̂)
    # Rearranging to match the answer choices format:
    # H_out = (K_0 / (1 + 2 * μ_0 / μ)) * (R^3/r^3) * (2*cos(θ)r̂ + sin(θ)θ̂)
    H_out_coeff = f"{K0} / (1 + (2*{mu0}/{mu}))"
    H_out_dipole_term = f"({R}^3/{r}^3) * (2*cos({theta}) {r_hat} + sin({theta}) {theta_hat})"

    print("The final expression for the magnetic field H(r, θ) is:")
    print("H(r, θ) = ")
    print("{")
    print(f"  {H_in_vector},  for 0 < r < R (inside)")
    print(f"  {H_out_coeff} * {H_out_dipole_term},  for r > R (outside)")
    print("}")
    print("\nThis matches answer choice E.")

solve_magnetic_field()