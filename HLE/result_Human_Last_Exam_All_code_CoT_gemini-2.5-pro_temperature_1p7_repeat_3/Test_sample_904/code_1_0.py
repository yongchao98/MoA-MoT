def solve_fluid_interface_equation():
    """
    This script derives the governing linear equation for the fluid interface
    and identifies the coefficients A(r) and B(r).
    """

    print("### Derivation of the Governing Equation ###\n")

    # Step 1: The Young-Laplace Equation
    print("Step 1: The Young-Laplace Equation")
    print("The relationship between pressure difference (delta_P) across a fluid interface,")
    print("surface tension (gamma), and mean curvature (H) is given by the Young-Laplace equation:")
    print("  delta_P = gamma * (2 * H)\n")

    # Step 2: Mean Curvature in Cylindrical Coordinates
    print("Step 2: Mean Curvature in Cylindrical Coordinates")
    print("For an axisymmetric surface defined by a function z = xi(r), the exact mean curvature H is:")
    print("  H = 0.5 * [ (d^2(xi)/dr^2) / (1 + (d(xi)/dr)^2)^(3/2) + (d(xi)/dr) / (r * (1 + (d(xi)/dr)^2)^(1/2)) ]\n")

    # Step 3: Linearization
    print("Step 3: Linearization for Small Displacements")
    print("The problem states that the displacement is small, which means the slope d(xi)/dr is also small.")
    print("Under this approximation, (d(xi)/dr)^2 can be neglected compared to 1.")
    print("The linearized mean curvature (H_lin) becomes:")
    print("  H_lin ≈ 0.5 * [ d^2(xi)/dr^2 + (1/r) * d(xi)/dr ]\n")

    # Step 4: Formulating the Pressure Balance
    print("Step 4: Formulating the Pressure Balance Equation")
    print("The system is in equilibrium, meaning the pressure exerted by the electric field (P_el) is balanced")
    print("by the capillary pressure from the surface tension. Gravity is negligible.")
    print("  P_el = delta_P_capillary = gamma * (2 * H_lin)")
    print("Substituting the linearized curvature, we get:")
    print("  P_el = gamma * [ d^2(xi)/dr^2 + (1/r) * d(xi)/dr ]\n")

    # Step 5: Rearranging to the Standard Form
    print("Step 5: Rearranging to the Standard Form")
    print("To match the requested form A(r)*xi'' + B(r)*xi' + C(r, xi) = 0, we rearrange the equation.")
    print("Let the electrostatic pressure term be represented by C(r, xi) = -P_el(r, xi).")
    print("The equation becomes:")
    print("  gamma * d^2(xi)/dr^2 + (gamma / r) * d(xi)/dr - P_el(r, xi) = 0\n")

    # Step 6: Identifying Coefficients A(r) and B(r)
    print("Step 6: Identifying Coefficients A(r) and B(r)")
    print("By comparing the derived equation with the standard form:")
    print("  A(r) * (d^2(xi)/dr^2) + B(r) * (d(xi)/dr) + C(r, xi) = 0")
    print("We can directly identify the coefficients A(r) and B(r).\n")

    # Final Result
    print("### Final Result ###")
    print("The coefficients of the governing linear equation are:")
    A_r = "gamma"
    B_r = "gamma / r"
    print(f"A(r) = {A_r}")
    print(f"B(r) = {B_r}")
    print("\nwhere 'gamma' is the surface tension between the two fluids and 'r' is the radial position.")

if __name__ == "__main__":
    solve_fluid_interface_equation()