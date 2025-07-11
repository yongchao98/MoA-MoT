def derive_interface_equation():
    """
    This script provides a step-by-step derivation of the coefficients A(r) and B(r)
    for the governing linear equation of a fluid interface in a cylindrical system.
    """

    print("--- Derivation of the Governing Equation for the Fluid Interface xi(r) ---\n")
    
    # --- Step 1: The Young-Laplace Equation ---
    print("Step 1: The Young-Laplace Equation")
    print("The shape of the interface is determined by the pressure balance across it.")
    print("The pressure jump due to surface tension (Laplace pressure, P_st) is given by:")
    print("  P_st = gamma * kappa")
    print("where 'gamma' is the surface tension and 'kappa' is the mean curvature of the interface.\n")

    # --- Step 2: Linearized Curvature in Cylindrical Coordinates ---
    print("Step 2: Linearized Curvature")
    print("For an interface with cylindrical symmetry described by z = xi(r), the mean curvature is:")
    print("  kappa = d^2(xi)/dr^2 * (1 + (d(xi)/dr)^2)^(-3/2) + (1/r) * d(xi)/dr * (1 + (d(xi)/dr)^2)^(-1/2)")
    print("Assuming small displacements, the slope d(xi)/dr is small. We can linearize this expression:")
    print("  (d(xi)/dr)^2 << 1")
    print("This simplifies the curvature to:")
    print("  kappa_linear = d^2(xi)/dr^2 + (1/r) * d(xi)/dr\n")

    # --- Step 3: Formulating the Pressure Balance Equation ---
    print("Step 3: Pressure Balance Equation")
    print("The system is in equilibrium when the Laplace pressure balances the electrostatic pressure, P_el(r).")
    print("Let's define the net pressure from external fields as C(r, xi). In this case, C is primarily due to the electric field.")
    print("The general force balance can be written as:")
    print("  P_laplace + C(r, xi) = 0")
    print("Substituting the linearized Laplace pressure gives:")
    print("  gamma * (d^2(xi)/dr^2 + (1/r) * d(xi)/dr) + C(r, xi) = 0\n")
    
    # --- Step 4: Identifying the Coefficients A(r) and B(r) ---
    print("Step 4: Identifying Coefficients")
    print("The derived equation is in the target form: A(r)*xi'' + B(r)*xi' + C(r, xi) = 0")
    print("  [gamma] * d^2(xi)/dr^2 + [gamma / r] * d(xi)/dr + C(r, xi) = 0")
    print("By comparing the terms, we can identify the coefficients A(r) and B(r).\n")

    # Define symbolic representations
    gamma_symbol = "\u03B3"  # Unicode for gamma
    r_symbol = "r"
    
    # Derived expressions for A(r) and B(r)
    A_r = f"{gamma_symbol}"
    B_r = f"{gamma_symbol} / {r_symbol}"

    print("--- Final Result ---")
    print(f"The coefficient A(r) is: {A_r}")
    print(f"The coefficient B(r) is: {B_r}\n")
    print(f"Here, '{gamma_symbol}' represents the constant surface tension and '{r_symbol}' is the radial coordinate.")

if __name__ == "__main__":
    derive_interface_equation()
