def derive_interface_equation():
    """
    This script explains the step-by-step derivation of the governing linear
    equation for the fluid interface shape ξ(r) and identifies the coefficients
    A(r) and B(r).
    """

    print("### Derivation of the Governing Equation for the Fluid Interface ###\n")
    print("The goal is to find the coefficients A(r) and B(r) in the equation:")
    print("A(r) * d²ξ/dr² + B(r) * dξ/dr + C(r, ξ) = 0\n")

    print("--- Step 1: The Young-Laplace Equation ---")
    print("The derivation starts with the Young-Laplace equation, which relates the pressure difference")
    print("across an interface (ΔP) to the surface tension (γ) and the mean curvature (H):")
    print("ΔP = γ * H")
    print("In this system, the pressure is due to the electric field, so ΔP = P_elec(r).\n")

    print("--- Step 2: Mean Curvature in Cylindrical Coordinates ---")
    print("For an axisymmetric surface described by z = ξ(r), the mean curvature is:")
    print("H = (1/r) * d/dr * [r * (dξ/dr) / sqrt(1 + (dξ/dr)²)]\n")

    print("--- Step 3: Linear Approximation for Small Displacements ---")
    print("Given that the interface displacement is small, the slope dξ/dr is also small.")
    print("This allows the approximation: sqrt(1 + (dξ/dr)²) ≈ 1.")
    print("The linearized curvature (H_lin) is therefore:")
    print("H_lin ≈ (1/r) * d/dr * [r * dξ/dr]\n")

    print("--- Step 4: Expanding the Linearized Curvature ---")
    print("Using the product rule for differentiation on the term [r * dξ/dr], we get:")
    print("d/dr[r * dξ/dr] = (1 * dξ/dr) + (r * d²ξ/dr²)")
    print("Substituting this back into the expression for H_lin:")
    print("H_lin ≈ (1/r) * [dξ/dr + r * d²ξ/dr²] = d²ξ/dr² + (1/r) * dξ/dr\n")

    print("--- Step 5: Assembling the Full Equation ---")
    print("Substituting H_lin and P_elec(r) into the Young-Laplace equation:")
    print("P_elec(r) = γ * (d²ξ/dr² + (1/r) * dξ/dr)")
    print("\nRearranging this to the standard form A(r)ξ'' + B(r)ξ' + C = 0 gives:")
    print("γ * d²ξ/dr² + (γ/r) * dξ/dr - P_elec(r) = 0\n")

    print("--- Step 6: Identifying the Coefficients ---")
    print("By comparing the derived equation to the target form, we can identify the coefficients.")
    
    # The coefficients are symbolic expressions
    a_r = "γ"
    b_r = "γ / r"
    
    print("\nThe coefficient A(r) multiplies the second derivative term (d²ξ/dr²):")
    print(f"A(r) = {a_r}")
    
    print("\nThe coefficient B(r) multiplies the first derivative term (dξ/dr):")
    print(f"B(r) = {b_r}")

# Run the derivation
if __name__ == "__main__":
    derive_interface_equation()
