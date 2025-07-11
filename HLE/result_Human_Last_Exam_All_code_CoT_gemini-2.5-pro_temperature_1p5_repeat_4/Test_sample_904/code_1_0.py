def derive_governing_equation():
    """
    This script derives the governing linear equation for the fluid interface ξ(r)
    and identifies the coefficients A(r) and B(r).
    """

    print("Derivation of the governing linear equation for the interfacial shape ξ(r):")
    print("=" * 75)

    # Step 1: Start with the Young-Laplace equation
    print("Step 1: The foundation of this problem is the Young-Laplace equation. It states that the pressure")
    print("difference (ΔP) across a fluid interface is proportional to its mean curvature (κ) and the")
    print("surface tension (γ):")
    print("  ΔP = γ * κ\n")

    # Step 2: Define the mean curvature in cylindrical coordinates and linearize
    print("Step 2: The interface is described by its height ξ as a function of the radial position r, i.e., ξ(r).")
    print("For an axisymmetric surface in cylindrical coordinates, the exact mean curvature is complex.")
    print("However, the problem states we can use a linear analysis (small displacements), which means the slope")
    print("dξ/dr is very small. In this limit, the mean curvature simplifies to:")
    print("  κ ≈ d²ξ/dr² + (1/r) * dξ/dr\n")

    # Step 3: Formulate the pressure balance
    print("Step 3: The pressure jump from surface tension (ΔP) must be balanced by any external pressures.")
    print("In this system, the external pressure is caused by the electric field, and gravity is negligible.")
    print("Let's denote this external pressure term as P_ext(r, ξ). The balance is:")
    print("  ΔP = P_ext(r, ξ)\n")

    # Step 4: Assemble the full equation
    print("Step 4: Now, we substitute the expressions for ΔP and κ into one equation. This gives the")
    print("relationship between the interface shape ξ(r) and the external pressure:")
    print("  γ * (d²ξ/dr² + (1/r) * dξ/dr) = P_ext(r, ξ)\n")

    # Step 5: Rearrange to the standard form
    print("Step 5: The final step is to arrange this equation into the standard linear form requested:")
    print("  A(r) * d²ξ/dr² + B(r) * dξ/dr + C(r, ξ) = 0\n")
    print("To do this, we simply move the external pressure term to the left side:")
    print("  γ * d²ξ/dr² + (γ/r) * dξ/dr - P_ext(r, ξ) = 0\n")

    # Step 6: Identify the coefficients A(r) and B(r)
    print("Step 6: By comparing our derived equation to the standard form, we can identify the coefficients.")
    print("The A(r) and B(r) terms arise from the surface tension effects, while the C(r, ξ) term")
    print("represents the influence of the external electric field (C(r, ξ) = -P_ext(r, ξ)).")
    print("-" * 75)
    
    # Define and print the final coefficients
    A_r = "γ"
    B_r = "γ / r"
    
    print("The identified coefficients are:")
    print(f"  A(r) = {A_r}")
    print(f"  B(r) = {B_r}")
    print("-" * 75)


if __name__ == "__main__":
    derive_governing_equation()