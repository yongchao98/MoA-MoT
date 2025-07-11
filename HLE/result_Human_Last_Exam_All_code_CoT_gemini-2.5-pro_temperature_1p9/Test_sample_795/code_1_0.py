import sympy

def solve_magnetization_curve():
    """
    This function derives and prints the analytical expression for the initial
    magnetization curve of the specified superconductor.
    """
    # Define the symbols for the physical quantities involved.
    # M: Magnetization
    # H: Applied magnetic field
    # a: Half-width of the superconductor bar
    # Jc: Critical current density
    M, H, a, Jc = sympy.symbols('M H a J_c')

    # The problem asks for the analytical expression for the initial magnetization curve.
    # Based on the Bean critical-state model for a slab of width 2a in a parallel field H,
    # the derivation steps are as follows:
    # 1. The field penetration depth 'd' is found to be d = H / Jc.
    # 2. The internal field profile H_y(x) is integrated across the slab width [-a, a] to find the average field <H_y>.
    #    The result is <H_y> = H^2 / (2 * a * Jc).
    # 3. Magnetization M is defined as M = <H_y> - H.
    # 4. Substituting the expression for <H_y>, we get the magnetization curve.

    magnetization_expression = H**2 / (2 * a * Jc) - H

    # The problem requests to output each number in the final equation.
    # We can write the expression as M = -1*H + 1*H**2 / (2*a*Jc)

    term1_coeff = -1
    term2_num_coeff = 1
    term2_den_coeff = 2

    print("The analytical expression for the initial magnetization M as a function of the applied field H is:")
    print(f"M(H) = ({term1_coeff})*H + ({term2_num_coeff})*H**2 / (({term2_den_coeff})*a*J_c)")

    # This expression is valid up to the full penetration field, H_p = a * Jc.
    print(f"\nThis expression is valid for an applied field H in the range: 0 <= H <= a*J_c.")

solve_magnetization_curve()

# The final analytical expression for the initial magnetization M(H) is:
# M(H) = -H * (1 - H / (2*a*Jc))
# Expanding this gives the form printed by the code.
final_answer_expr = "-H * (1 - H / (2*a*Jc))"
<<<M(H) = -H + H**2/(2*a*Jc)>>>