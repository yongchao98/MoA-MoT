import sympy

def solve_ghz_purification():
    """
    This script analytically solves for the product of the output fidelity 
    and success probability in a GHZ state purification protocol.
    It uses the sympy library to perform symbolic calculations and present the final formula.
    """
    
    # --- Introduction and Plan ---
    print("This script calculates the product of the output fidelity (F_out) and the success probability (P_succ) for a given GHZ state purification protocol.")
    print("The problem is solved analytically by deriving a formula in terms of the input fidelities F1 and F2.\n")
    print("Plan:")
    print("1. Express the total input state as a linear combination of four basic operator terms.")
    print("2. Calculate the contribution of each term to the final quantity (F_out * P_succ) after the protocol.")
    print("3. Sum the contributions, weighted by their coefficients, to get the final expression.")
    print("4. Simplify the expression to get the final formula.\n")

    # --- Symbolic Calculation ---

    # Step 1: Define symbolic variables and coefficients
    print("--- Step 1: Define symbolic variables and coefficients ---")
    F1, F2 = sympy.symbols('F1 F2')
    
    # Coefficients from the state definitions
    # rho_GHZ(F1) = alpha * |GHZ><GHZ| + beta * I_3
    # rho_Bell(F2) = gamma * |Phi+><Phi+| + delta * I_2
    alpha = (8*F1 - 1) / 7
    beta = (1 - F1) / 7
    gamma = (4*F2 - 1) / 3
    delta = (1 - F2) / 3
    
    print("The input state rho_in can be written as:")
    print("rho_in = (a * |GHZ><GHZ| + b * I_3) ⊗ (c * |Phi+><Phi+| + d * I_2), where:")
    print(f"a = (8*F1 - 1)/7")
    print(f"b = (1 - F1)/7")
    print(f"c = (4*F2 - 1)/3")
    print(f"d = (1 - F2)/3\n")
    
    # Step 2: Calculate contributions for each of the four terms
    print("--- Step 2: Calculate contributions for the four basis terms ---")
    print("The desired quantity, F_out * P_succ, is linear. We calculate its value (V_i) for each of the four basic terms resulting from the expansion of rho_in:")
    # These values (V1, V2, V3, V4) are derived from the physics of the protocol as explained in the thought process.
    # V_i = Tr(Pi_out * U * rho_term_i * U_dagger)
    V1 = 1  # For |GHZ><GHZ| ⊗ |Phi+><Phi+|
    V2 = 1  # For |GHZ><GHZ| ⊗ I_2
    V3 = 1  # For I_3 ⊗ |Phi+><Phi+|
    V4 = 2  # For I_3 ⊗ I_2
    print(f"V1 (for |GHZ><GHZ| ⊗ |Phi+><Phi+|): {V1}")
    print(f"V2 (for |GHZ><GHZ| ⊗ I_2):         {V2}")
    print(f"V3 (for I_3 ⊗ |Phi+><Phi+|):        {V3}")
    print(f"V4 (for I_3 ⊗ I_2):                {V4}\n")

    # Step 3 & 4: Combine the results and simplify
    print("--- Step 3 & 4: Combine contributions and simplify ---")
    # Total expression is the sum of contributions weighted by coefficients
    # F_out * P_succ = alpha*gamma*V1 + alpha*delta*V2 + beta*gamma*V3 + beta*delta*V4
    total_expr = alpha * gamma * V1 + alpha * delta * V2 + beta * gamma * V3 + beta * delta * V4
    
    print("The full expression for F_out * P_succ is:")
    print("a*c*V1 + a*d*V2 + b*c*V3 + b*d*V4\n")
    
    # Simplify the final expression
    simplified_expr = sympy.simplify(total_expr)
    
    # Extract numerator and denominator to print clearly
    num, den = sympy.fraction(simplified_expr)
    expanded_num = sympy.expand(num)
    
    # Extract coefficients to fulfill the "output each number" requirement
    term_f1f2 = expanded_num.coeff(F1 * F2)
    # The .coeff(var, 1) method can sometimes be tricky if a variable is not present.
    # This safer approach extracts the coefficient of F1 from the part of the expression not containing F2.
    term_f1 = sympy.expand(expanded_num.coeff(F2, 0)).coeff(F1, 1)
    term_f2 = sympy.expand(expanded_num.coeff(F1, 0)).coeff(F2, 1)
    const_term = sympy.expand(expanded_num.coeff(F1,0).coeff(F2,0))

    print("The final simplified equation for the product of the output fidelity and success probability is:")
    # Print the equation with each number explicitly shown
    final_equation_str = f"({term_f1f2}*F1*F2 + ({term_f1})*F1 + ({term_f2})*F2 + {const_term}) / {den}"
    print(final_equation_str)

if __name__ == "__main__":
    solve_ghz_purification()
