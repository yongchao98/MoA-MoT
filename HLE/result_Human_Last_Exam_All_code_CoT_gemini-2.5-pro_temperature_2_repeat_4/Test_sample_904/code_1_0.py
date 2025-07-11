def get_governing_equation_coefficients():
    """
    This script presents the derived coefficients A(r) and B(r) for the
    governing linear equation of the fluid interface shape ξ(r).

    The equation is of the form:
    A(r) * d²ξ/dr² + B(r) * dξ/dr + C(r, ξ) = 0
    """

    # Symbolic representations for physical quantities
    gamma = "γ"  # Surface tension
    r = "r"      # Radial position

    # The derived symbolic expressions for the coefficients
    A_r_expression = gamma
    B_r_expression = f"{gamma} / {r}"

    # Print the final derived expressions clearly
    print("For the linearized governing equation of the interfacial shape ξ(r):")
    print("A(r) * d²ξ/dr² + B(r) * dξ/dr + C(r, ξ) = 0")
    print("-" * 50)
    
    print("The coefficient A(r) is:")
    print(f"A(r) = {A_r_expression}")
    
    print("\nThe coefficient B(r) is:")
    print(f"B(r) = {B_r_expression}")
    
    print("-" * 50)
    print("Note: These coefficients are derived from the linearized surface tension term (Laplace pressure) in cylindrical coordinates.")

if __name__ == "__main__":
    get_governing_equation_coefficients()
