import sympy

def solve_purification_protocol():
    """
    This function calculates the product of the successful output fidelity 
    and the success probability for the given GHZ state purification protocol.
    The calculation is performed symbolically in terms of the input fidelities F1 and F2.
    """
    
    # Define symbolic variables for the fidelities F1 and F2
    F1, F2 = sympy.symbols('F1 F2')

    # Define the coefficients from the input state definitions as per the problem statement
    # For the 3-qubit GHZ state: rho_GHZ(F1) = alpha1 * |GHZ><GHZ| + beta1 * I_8
    alpha1 = (8*F1 - 1) / 7
    beta1 = (1 - F1) / 7
    
    # For the 2-qubit Bell state: rho_Bell(F2) = alpha2 * |Phi+><Phi+| + beta2 * I_4
    alpha2 = (4*F2 - 1) / 3
    beta2 = (1 - F2) / 3

    # As derived from the analysis of the quantum circuit, the target quantity 
    # <GHZ|tilde_rho_succ|GHZ> is calculated for each of the four components
    # of the expanded input density matrix.
    # Q1 is for the component: |GHZ><GHZ| (x) |Phi+><Phi+|
    # Q2 is for the component: |GHZ><GHZ| (x) I_4
    # Q3 is for the component: I_8 (x) |Phi+><Phi+|
    # Q4 is for the component: I_8 (x) I_4
    Q1 = 1
    Q2 = 1
    Q3 = 1
    Q4 = 2

    # The total product is the linear combination of these quantities,
    # weighted by their respective coefficients from the density matrix expansion.
    total_product = alpha1 * alpha2 * Q1 + alpha1 * beta2 * Q2 + beta1 * alpha2 * Q3 + beta1 * beta2 * Q4

    # Use sympy to simplify the final expression
    simplified_product = sympy.simplify(total_product)

    # To display the final formula with each number explicitly, we extract the
    # numerator and denominator of the resulting rational function.
    num, den = simplified_product.as_numer_denom()

    # We then create a polynomial from the numerator to extract the coefficients.
    num_poly = sympy.Poly(num, F1, F2)
    c_f1f2 = num_poly.coeff_monomial(F1*F2)
    c_f1 = num_poly.coeff_monomial(F1)
    c_f2 = num_poly.coeff_monomial(F2)
    c_const = num_poly.coeff_monomial(1)
    denominator = den

    print("The product of the successful output fidelity and the success probability is a function of the input fidelities F1 and F2.")
    print("The final derived expression is:")
    print(f"({c_f1f2}*F1*F2 + ({c_f1})*F1 + ({c_f2})*F2 + {c_const}) / {denominator}")

if __name__ == '__main__':
    solve_purification_protocol()