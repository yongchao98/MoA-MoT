import sympy

def calculate_energy_shift():
    """
    Calculates the second-order energy shift for a hydrogen atom
    due to the relativistic kinetic energy correction.
    """
    # Define the principal and angular momentum quantum numbers for the state
    n_val = 3
    l_val = 2
    # Define the atomic number for Hydrogen
    Z_val = 1

    # Use sympy for exact fractional arithmetic
    n, l, Z = sympy.Integer(n_val), sympy.Integer(l_val), sympy.Integer(Z_val)

    # The formula for the second-order energy shift due to the relativistic kinetic
    # energy correction (H' = -p^4 / (8*m^3*c^2)) can be found in advanced
    # quantum mechanics texts. For l >= 1, a very good approximation is:
    # ΔE = - (m*c^2 * (Z*α)^6 / n^3) * [ 1/(2ℓ+1)^2 - 11/(32n(2ℓ+1)) + 3/(4n^2(2ℓ+1)^2) ]
    
    # We will calculate the final coefficient step-by-step.
    
    print("The formula for the second-order energy shift is:")
    print("ΔE = - (m*c^2 * (Z*α)**6 / n**3) * [A + B + C]")
    print(f"\nFor a hydrogen atom (Z={Z_val}) in the state n={n_val}, ℓ={l_val}:")
    
    # --- Part 1: Calculate the terms inside the square brackets ---
    print("\nStep 1: Calculate the terms A, B, and C in the square brackets.")
    
    # Term A: 1/(2ℓ+1)^2
    term_A_num = 1
    term_A_den = (2*l + 1)**2
    term_A = sympy.Rational(term_A_num, term_A_den)
    print(f"Term A = 1 / (2*{l_val} + 1)**2 = {term_A_num} / {term_A_den}")
    
    # Term B: -11/(32n(2ℓ+1))
    term_B_num = -11
    term_B_den = 32 * n * (2*l + 1)
    term_B = sympy.Rational(term_B_num, term_B_den)
    print(f"Term B = -11 / (32 * {n_val} * (2*{l_val} + 1)) = {term_B_num} / {term_B_den}")

    # Term C: 3/(4n^2(2ℓ+1)^2)
    term_C_num = 3
    term_C_den = 4 * n**2 * (2*l + 1)**2
    term_C = sympy.Rational(term_C_num, term_C_den)
    print(f"Term C = 3 / (4 * {n_val}**2 * (2*{l_val} + 1)**2) = {term_C_num} / {term_C_den}")

    # Sum of the terms in the bracket
    bracket_term = term_A + term_B + term_C
    print(f"\nSum of bracket terms [A + B + C] = {term_A} + ({term_B}) + {term_C} = {bracket_term}")

    # --- Part 2: Calculate the prefactor ---
    print("\nStep 2: Calculate the prefactor outside the brackets.")
    prefactor_coeff_num = -(Z**6)
    prefactor_coeff_den = n**3
    prefactor_coeff = sympy.Rational(prefactor_coeff_num, prefactor_coeff_den)
    print(f"Prefactor coefficient = -(Z**6 / n**3) = -({Z_val}**6 / {n_val}**3) = {prefactor_coeff}")
    
    # --- Part 3: Combine to get the final result ---
    print("\nStep 3: Multiply the prefactor coefficient by the bracket sum.")
    total_coeff = prefactor_coeff * bracket_term
    print(f"Total coefficient = ({prefactor_coeff}) * ({bracket_term}) = {total_coeff}")
    
    # --- Final Answer ---
    print("\nThe final expression for the energy shift, ΔE, is:")
    # Using unicode characters for a nicer print format
    m_sym, c_sym, alpha_sym = sympy.symbols('m c α')
    final_expression = total_coeff * m_sym * c_sym**2 * alpha_sym**6
    
    # Print the equation with all the numbers
    print(f"ΔE = ({total_coeff.p}/{total_coeff.q}) * m*c**2 * α**6")

# Execute the calculation
calculate_energy_shift()