import sympy

def calculate_relativistic_correction():
    """
    Calculates the second-order energy shift for a hydrogen atom
    in the n=3, l=2 state due to the relativistic kinetic energy correction.
    """
    # Define symbols for fundamental constants and quantum numbers
    m, c, alpha = sympy.symbols('m c alpha', positive=True, real=True)
    n_val, l_val = 3, 2

    # Rydberg energy in terms of fundamental constants
    E_R = sympy.Rational(1, 2) * m * c**2 * alpha**2

    # --- Step 1: Calculate the first-order energy correction E^(1) ---

    # Bohr energy for the given n
    E_n0 = -E_R / n_val**2
    
    # Calculate the term in the parenthesis for E^(1)
    factor1 = (4 * n_val) / (l_val + sympy.Rational(1, 2)) - 3
    
    # Calculate E^(1)
    E_1 = - (E_n0**2) / (2 * m * c**2) * factor1
    
    # Simplify the expression for E^(1) in terms of E_R and alpha
    E_1_simplified = sympy.simplify(E_1.subs(m * c**2, 2 * E_R / alpha**2))
    
    print(f"Calculation for the first-order shift E^(1) for n={n_val}, l={l_val}:")
    print(f"E_1 = - ((-E_R/{n_val**2})^2) / (2*m*c^2) * ( (4*{n_val})/({l_val}+1/2) - 3 )")
    print(f"E_1 = {E_1_simplified}")
    print("-" * 30)

    # --- Step 2: Calculate the second-order energy correction E^(2) ---

    # Use the formula relating E^(2) to E^(1)
    E_2 = -(E_1**2) * n_val**2 / (E_R * (2 * l_val + 1))
    
    # Simplify the final expression
    # Substitute the simplified E_1 first
    E_2_intermediate = -(E_1_simplified**2) * n_val**2 / (E_R * (2 * l_val + 1))
    E_2_final_rydberg = sympy.simplify(E_2_intermediate)
    
    # Express the final result in terms of m, c, and alpha
    E_2_final_fundamental = sympy.simplify(E_2_final_rydberg.subs(E_R, sympy.Rational(1, 2) * m * c**2 * alpha**2))

    print(f"Calculation for the second-order shift E^(2):")
    # To show the equation clearly, we use the simplified E_1 expression
    e1_sq_term = f"({E_1_simplified})**2"
    denom_term = f"(E_R * (2*{l_val}+1))"
    num_term = f"{n_val**2}"
    print(f"E_2 = - {e1_sq_term} * ( {num_term} / {denom_term} )")
    print(f"E_2 = {E_2_final_rydberg}")
    print(f"E_2 = {E_2_final_fundamental}")

    # Extract the numerical coefficient for the final answer
    coefficient = E_2_final_fundamental / (m * c**2 * alpha**6)
    print("\nFinal Answer:")
    print(f"The second-order energy shift is E^(2) = ({coefficient}) * m*c^2*alpha^6")
    
# Execute the calculation
calculate_relativistic_correction()