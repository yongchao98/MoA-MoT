import sympy

def find_fixed_point_coupling():
    """
    Calculates and prints the leading order expression for the Wilson-Fisher
    fixed point coupling u* in phi^4 theory near d=4 dimensions.
    """
    # --- Introduction and Setup ---
    print("This script finds the fixed point coupling u^* for the phi^4 theory in d = 4 - epsilon dimensions.")
    print("The calculation uses the one-loop renormalization group (RG) beta function.")
    print("\nThe beta function, beta(u), describes the running of the dimensionless coupling 'u' with energy scale.")
    print("To leading order in epsilon and u, it is given by:")
    print("beta(u) = -epsilon * u + B * u^2")

    # Define the coefficient B for single-component phi^4 theory
    # This corresponds to the interaction Lagrangian L_int = (lambda/4!) * phi^4
    # The value of B is a standard result from one-loop calculations.
    B_num = 3
    B_den_coeff = 16
    
    print(f"\nFor the standard single-component phi^4 theory, the constant B = {B_num} / ({B_den_coeff} * pi^2).")

    print("\nA fixed point u^* is a value of the coupling where the theory is scale-invariant,")
    print("which means the beta function is zero: beta(u^*) = 0.")
    print("\nSolving the equation: -epsilon * u^* + B * (u^*)^2 = 0")
    print("We are looking for the non-trivial solution where u^* is not zero.")

    # --- Symbolic Calculation ---
    
    # Define the symbols needed for the equation
    u_star = sympy.Symbol('u^*')
    epsilon = sympy.Symbol('epsilon')
    pi = sympy.pi

    # Define the constant B
    B = B_num / (B_den_coeff * pi**2)

    # Define the beta function equation
    beta_function_eq = sympy.Eq(-epsilon * u_star + B * u_star**2, 0)
    
    # Solve the equation for u*
    solutions = sympy.solve(beta_function_eq, u_star)

    # The solutions are [0, epsilon/B]. We want the non-trivial one.
    wilson_fisher_fixed_point = None
    for sol in solutions:
        if sol != 0:
            wilson_fisher_fixed_point = sol
            break
            
    # --- Print the Result ---
    print("\n-------------------------------------------")
    print("The leading order expression for the fixed point coupling is:")
    
    # The solution is epsilon / B = epsilon / (3 / (16*pi^2)) = (16*pi^2 / 3) * epsilon
    numerator_coeff = B_den_coeff
    denominator = B_num
    
    # Final expression in a formatted string showing all numbers
    final_expression_str = f"{u_star.name} = ({numerator_coeff} * {pi.name}^2 / {denominator}) * {epsilon.name}"
    
    print(final_expression_str)
    print("-------------------------------------------")

find_fixed_point_coupling()