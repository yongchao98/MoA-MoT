import sympy

def solve_purification_product():
    """
    This function calculates and prints the product of the success probability and
    the output fidelity for the described GHZ state purification protocol.
    The calculation is performed symbolically using the sympy library.
    """
    # Step 1: Explain the methodology.
    print("Methodology:")
    print("The product of the output fidelity (F_out) and success probability (P_succ) is calculated.")
    print("This product, F_out * P_succ, is equivalent to Tr(|GHZ><GHZ| * rho_out_unnormalized),")
    print("where rho_out_unnormalized is the unnormalized output state on the first 3 qubits.")
    print("The input state is decomposed into four terms, and the contribution from each is calculated and summed.")
    print("-" * 20)

    # Step 2: Define symbolic variables for the fidelities F1 and F2.
    F1, F2 = sympy.symbols('F1 F2')

    # Step 3: Define the coefficients p1, q1, p2, q2 based on the problem description.
    # rho_GHZ(F1) = p1*|GHZ><GHZ| + q1*I_8
    # rho_Bell(F2) = p2*|Phi+><Phi+| + q2*I_4
    p1 = (8 * F1 - 1) / 7
    q1 = (1 - F1) / 7
    p2 = (4 * F2 - 1) / 3
    q2 = (1 - F2) / 3

    # Step 4: The value of Tr(|GHZ><GHZ| * rho_out_unnormalized) has been derived for each of the four components:
    # 1. For |GHZ><GHZ| x |Bell><Bell| component (coeff p1*p2), the value is 1.
    # 2. For |GHZ><GHZ| x I_2 component (coeff p1*q2), the value is 1.
    # 3. For I_3 x |Bell><Bell| component (coeff q1*p2), the value is 1.
    # 4. For I_3 x I_2 component (coeff q1*q2), the value is 2.
    # The total product is the sum of these contributions, weighted by their coefficients.
    product_expr = p1 * p2 * 1 + p1 * q2 * 1 + q1 * p2 * 1 + q1 * q2 * 2

    # Step 5: Simplify the expression using sympy.
    simplified_expr = sympy.simplify(product_expr)

    # Step 6: Format and print the final equation with each number explicitly shown.
    print("The product of the output fidelity and the success probability is:")
    
    # Get the numerator and denominator of the simplified expression
    num, den = simplified_expr.as_numer_den()
    
    # Represent the numerator as a polynomial in F1 and F2 to extract coefficients
    poly_num = sympy.Poly(num, F1, F2)
    
    terms = []
    # Iterate through the terms of the polynomial (coefficient, powers of variables)
    for monom, coeff in poly_num.terms():
        (p_F1, p_F2) = monom
        
        term_str = ""
        # Handle the coefficient
        if coeff == 1 and (p_F1 > 0 or p_F2 > 0):
            term_str += "" # No '1' for terms like F1*F2
        elif coeff == -1 and (p_F1 > 0 or p_F2 > 0):
            term_str += "-"
        else:
            term_str += str(coeff)
            if p_F1 > 0 or p_F2 > 0:
                 term_str += "*"

        # Handle the variables
        if p_F1 > 0:
            term_str += "F1"
            if p_F1 > 1:
                term_str += f"**{p_F1}"
            if p_F2 > 0:
                term_str += "*"
        if p_F2 > 0:
            term_str += "F2"
            if p_F2 > 1:
                term_str += f"**{p_F2}"
        
        terms.append(term_str)

    # Join terms with ' + ' and clean up for readability (e.g., ' + -' -> ' - ')
    equation_num_str = " + ".join(terms).replace(" + -", " - ")

    final_equation = f"({equation_num_str}) / {den}"
    print(final_equation)
    
    # Final answer in the required format
    print(f"\n<<< {simplified_expr} >>>")

# Execute the function
solve_purification_product()