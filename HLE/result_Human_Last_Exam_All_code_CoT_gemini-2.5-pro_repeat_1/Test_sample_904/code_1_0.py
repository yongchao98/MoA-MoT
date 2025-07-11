def solve_fluid_equation():
    """
    This script derives the governing linear equation for the interfacial shape xi(r)
    and identifies the coefficients A(r) and B(r).
    """

    # Explanation of symbols used in the derivation
    # gamma: Surface tension between the two fluids
    # r: Radial coordinate
    # xi(r): Displacement of the fluid interface at radius r
    # xi'(r): First derivative of xi with respect to r (d(xi)/dr)
    # xi''(r): Second derivative of xi with respect to r (d^2(xi)/dr^2)
    # kappa: Mean curvature of the interface
    # Delta_P: Pressure difference across the interface due to the electric field

    print("Derivation of the Governing Linear Equation for the Fluid Interface\n")

    print("Step 1: The Young-Laplace Equation")
    print("The shape of the interface is governed by the Young-Laplace equation, which balances the pressure difference (Delta_P) across the interface with the capillary pressure due to surface tension (gamma) and curvature (kappa).")
    print("Equation: Delta_P = gamma * kappa\n")

    print("Step 2: Curvature in Cylindrical Coordinates")
    print("For an axisymmetric surface described by z = xi(r), the exact mean curvature kappa is given by:")
    print("kappa = xi''(r) / (1 + xi'(r)^2)^(3/2)  +  xi'(r) / (r * (1 + xi'(r)^2)^(1/2))\n")

    print("Step 3: Linearization for Small Displacements")
    print("The problem states that a linear analysis is valid, which implies the slope of the interface is small (xi'(r) << 1).")
    print("Under this approximation, (1 + xi'(r)^2) ≈ 1. The curvature expression simplifies significantly:")
    print("kappa_linear ≈ xi''(r) + (1/r) * xi'(r)\n")

    print("Step 4: Formulating the Linear Differential Equation")
    print("Substituting the linearized curvature back into the Young-Laplace equation gives:")
    print("Delta_P = gamma * (xi''(r) + (1/r) * xi'(r))\n")

    print("Step 5: Rearranging to the Target Form")
    print("The problem asks for the equation in the form: A(r)*xi''(r) + B(r)*xi'(r) + C(r, xi(r)) = 0.")
    print("We can rearrange our equation by moving all terms to one side:")
    print("gamma * xi''(r) + (gamma / r) * xi'(r) - Delta_P = 0\n")

    print("Step 6: Identifying the Coefficients A(r) and B(r)")
    print("By comparing our derived equation with the target form, we can directly identify the coefficients A(r) and B(r).")
    print("The term C(r, xi(r)) corresponds to the negative of the pressure difference, -Delta_P, which includes all effects from the electric field.\n")

    # Final Answer
    A_r = "gamma"
    B_r = "gamma / r"

    print("--- FINAL RESULT ---")
    print("The coefficients of the governing linear equation are:")
    print(f"A(r) = {A_r}")
    print(f"B(r) = {B_r}")

if __name__ == '__main__':
    solve_fluid_equation()