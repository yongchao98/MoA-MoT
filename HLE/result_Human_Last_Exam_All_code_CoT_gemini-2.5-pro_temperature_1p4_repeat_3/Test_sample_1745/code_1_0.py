import sympy as sp

def solve_edl_potential():
    """
    This function derives and prints the symbolic expression for the
    Electrical Double-Layer (EDL) potential distribution psi(y) based on the
    provided problem description.
    """
    # Define the symbols used in the equation
    y, z_1, beta, k, H = sp.symbols('y z_1 beta k H')
    psi = sp.Function('psi')(y)

    # According to the plan, the final derived expression for the potential distribution is:
    # psi(y) = z_1*(1 + beta*k) * sinh(k*(H/2 - y)) / sinh(k*H)
    # We will construct this expression to print it out clearly, ensuring
    # every variable and constant from the formula is shown.

    # Part 1: The zeta potential at the bottom wall
    term1 = f"z_1*(1 + beta*k)"

    # Part 2: The hyperbolic sine term in the numerator
    term2 = f"sinh(k*(H/2 - y))"

    # Part 3: The hyperbolic sine term in the denominator
    term3 = f"sinh(k*H)"

    # Print the final equation in a clear, single-line format
    print("The final expression for the Electrical double-layer potential distribution is:")
    # The following line prints the complete symbolic equation as requested.
    print(f"psi(y) = ({term1}) * {term2} / {term3}")

if __name__ == '__main__':
    solve_edl_potential()