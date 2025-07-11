import sympy as sp

def calculate_inductance_change():
    """
    This script calculates and displays the expression for the change in mutual
    inductance per unit length between the two circuits when the concentrator
    is added, in the limit where d >> h.
    """

    # Define the symbolic variables for the physical quantities.
    # mu_0: Permeability of free space
    # h: Separation of wires within a circuit
    # d: Separation between the centers of the two circuits
    # R1: Inner radius of the concentrator shell
    mu_0, h, d, R1 = sp.symbols('mu_0 h d R_1', real=True, positive=True)

    # Based on the physical analysis using the method of images and dipole
    # approximation, the change in mutual inductance per unit length (Delta M')
    # is derived. We construct the expression below.

    # The final expression is of the form:
    #   -8 * mu_0 * h**2 * R1**2
    # -----------------------------
    #   pi * (d**2 + 4 * R1**2)**2

    # Let's build and print this expression.
    
    # Numerical constants from the derivation
    const_num = -8
    const_den = 4
    exponent = 2

    # Build the numerator and denominator
    numerator = const_num * mu_0 * h**exponent * R1**exponent
    denominator = sp.pi * (d**exponent + const_den * R1**exponent)**exponent

    # The full expression
    delta_M_per_length = numerator / denominator

    # --- Output the results ---
    print("The expression for the change in mutual inductance per unit length (Delta M' = M2' - M1') is:")
    sp.pprint(delta_M_per_length, use_unicode=True)

    print("\nTo clarify the numbers in the final equation, it can be written as:")
    final_expression_str = f"({const_num}) * mu_0 * h**{exponent} * R1**{exponent} / (pi * (d**{exponent} + ({const_den}) * R1**{exponent})**{exponent})"
    print(final_expression_str)

if __name__ == '__main__':
    calculate_inductance_change()