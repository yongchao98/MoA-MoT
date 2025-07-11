def solve_magnetic_field():
    """
    This function outlines the derivation of the magnetic field for the given spherical shell problem
    and prints the final result in a step-by-step formatted way.
    """

    print("### Derivation of the Magnetic Field ###")
    print("\n--- Step 1: Define Magnetic Scalar Potential ---")
    print("In regions with no current, H = -∇Φ_m, where ∇²Φ_m = 0.")
    print("Due to the current's form K = K₀sin(θ)φ_hat, the solution involves the l=1 Legendre polynomial, P₁(cosθ) = cosθ.")

    print("\n--- Step 2: General Solutions for Potential ---")
    print("Inside (r < R), the potential must be finite at r=0:")
    print("  Φ_in(r, θ) = A₁ * r * cos(θ)")
    print("Outside (r > R), the potential must vanish as r → ∞:")
    print("  Φ_out(r, θ) = B₁ * r⁻² * cos(θ)")

    print("\n--- Step 3: Field Components from Potentials ---")
    print("H = -∇Φ_m = -(∂Φ/∂r r_hat + (1/r)∂Φ/∂θ θ_hat)")
    print("  H_in = -A₁cos(θ) r_hat + A₁sin(θ) θ_hat")
    print("  H_out = 2B₁r⁻³cos(θ) r_hat + B₁r⁻³sin(θ) θ_hat")

    print("\n--- Step 4: Apply Boundary Conditions at r=R ---")
    print("1. B_normal is continuous: μ*H_in,r = μ₀*H_out,r")
    print("   μ * (-A₁) = μ₀ * (2B₁R⁻³)  ==>  A₁ = -(2μ₀/μ) * B₁R⁻³   (Eq. 1)")
    print("2. H_tangential is discontinuous: H_out,θ - H_in,θ = K_φ")
    print("   B₁R⁻³ - A₁ = K₀           (Eq. 2)")

    print("\n--- Step 5: Solve for Coefficients A₁ and B₁ ---")
    print("Substitute (Eq. 1) into (Eq. 2):")
    print("  B₁R⁻³ - (-(2μ₀/μ) * B₁R⁻³) = K₀")
    print("  B₁R⁻³ * (1 + 2μ₀/μ) = K₀  ==> B₁ = (K₀*R³) / (1 + 2μ₀/μ)")
    print("And subsequently:")
    print("  A₁ = -(2μ₀/μ) * K₀ / (1 + 2μ₀/μ)")

    print("\n--- Step 6: Construct the Final Magnetic Field Expressions ---")
    print("Substitute A₁ and B₁ back into the field component equations.")

    # Construct the inside field expression
    print("\nFor the region inside (0 < r < R):")
    # H_in = -A₁ * (cos(θ)r_hat - sin(θ)θ_hat) = -A₁ * z_hat
    # H_in = - (-(2μ₀/μ) * K₀ / (1 + 2μ₀/μ)) * z_hat
    inside_field_coeff_numerator = "2 * μ₀"
    inside_field_coeff_denominator = f"μ * (1 + (2 * μ₀ / μ))"
    inside_field_coeff = f"({inside_field_coeff_numerator} / {inside_field_coeff_denominator}) * K₀"
    # Alternative form used in the options
    inside_field_coeff_alt = f"(2 * μ₀ / μ) * (K₀ / (1 + (2 * μ₀ / μ)))"

    print(f"H_inside(r, θ) = {inside_field_coeff_alt} * z_hat")
    print("\nNote that z_hat = cos(θ)r_hat - sin(θ)θ_hat.")
    print("This corresponds to a uniform field in the +z direction.")

    # Construct the outside field expression
    print("\nFor the region outside (R < r < ∞):")
    outside_field_coeff = f"K₀ / (1 + (2 * μ₀ / μ))"
    outside_field_radial = f"(R³/r³)"
    outside_field_angular = f"(2*cos(θ) r_hat + sin(θ) θ_hat)"
    print(f"H_outside(r, θ) = ({outside_field_coeff}) * {outside_field_radial} * {outside_field_angular}")
    print("\nThis is the field of a magnetic dipole.")

    print("\n--- Conclusion ---")
    print("Comparing these derived results with the provided choices, the correct one is Option E.")

solve_magnetic_field()
print("<<<E>>>")