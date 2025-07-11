def solve_potential_distribution():
    """
    This function formulates and prints the final expression for the
    Electrical double-layer (EDL) potential distribution, psi(y).
    
    The derivation is based on solving the linearized Poisson-Boltzmann equation
    d^2(psi)/dy^2 = k^2 * psi
    with the following boundary conditions:
    1. At y = 0: psi(0) = z_a1 = z_1 * (1 + beta * k)
    2. At y = H: psi(H) = z_a2 = 0 (since z_2 = 0)
    
    The final expression is printed as a formatted string.
    """

    # Define the symbols as strings for the final expression
    psi_y = "psi(y)"
    z_1 = "z_1"
    one = "1"
    beta = "beta"
    k = "k"
    H = "H"
    y = "y"

    # Construct the components of the expression
    coefficient = f"{z_1}*({one} + {beta}*{k})"
    sinh_numerator = f"sinh({k}*({H} - {y}))"
    sinh_denominator = f"sinh({k}*{H})"

    # Combine all parts into the final expression string
    final_expression = f"{psi_y} = {coefficient} * ({sinh_numerator} / {sinh_denominator})"
    
    print("The final expression for the Electrical double-layer potential distribution psi(y) is:")
    print(final_expression)

solve_potential_distribution()
<<<psi(y) = z_1*(1 + beta*k) * (sinh(k*(H - y)) / sinh(k*H))>>>