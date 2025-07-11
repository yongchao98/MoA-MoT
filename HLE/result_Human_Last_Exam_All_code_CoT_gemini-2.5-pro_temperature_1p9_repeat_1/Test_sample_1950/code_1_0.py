import sympy

def solve_purification_protocol():
    """
    Calculates and prints the product of success probability and output fidelity
    for the described GHZ state purification protocol.
    """
    # Define symbolic variables for the fidelities
    F1, F2 = sympy.symbols('F1 F2')

    # Define the coefficients from the input state definitions
    # For rho_GHZ(F1)
    alpha1 = (8*F1 - 1) / 7
    beta1 = (1 - F1) / 7
    # For rho_Bell(F2)
    alpha2 = (4*F2 - 1) / 3
    beta2 = (1 - F2) / 3

    # The contributions (S-values) calculated for the four basic operator products
    # S1: for |GHZ><GHZ| (x) |Phi+><Phi+|
    S1 = sympy.Rational(1, 4)
    # S2: for |GHZ><GHZ| (x) I_4
    S2 = 1
    # S3: for I_8 (x) |Phi+><Phi+|
    S3 = 1
    # S4: for I_8 (x) I_4
    S4 = 2

    # Calculate the total product S by summing the weighted contributions
    total_product = (alpha1 * alpha2 * S1 +
                     alpha1 * beta2 * S2 +
                     beta1 * alpha2 * S3 +
                     beta1 * beta2 * S4)

    # Simplify the final symbolic expression
    simplified_product = sympy.simplify(total_product)
    
    # Extract numerator and denominator for formatted printing
    num, den = sympy.fraction(simplified_product)
    
    # Extract coefficients from the numerator polynomial
    # Note: .coeff() method returns 0 if the term doesn't exist.
    poly_num = sympy.Poly(num, F1, F2)
    c_f1 = poly_num.coeff_monomial(F1)
    c_f2 = poly_num.coeff_monomial(F2)
    c_f1f2 = poly_num.coeff_monomial(F1*F2)
    # The constant term is the one with F1=0, F2=0
    c_const = poly_num.subs({F1: 0, F2: 0})

    # Print the explanation and the final equation
    print("The product of the successful output fidelity and the success probability (S) is given by the formula:")
    print(f"S = ({c_f1}*F1 + {c_f2}*F2 + {c_f1f2}*F1*F2 + {c_const}) / {den}")
    # Also printing the raw sympy result for verification
    # print("\nRaw symbolic expression:")
    # print(simplified_product)


if __name__ == '__main__':
    solve_purification_protocol()