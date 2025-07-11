import sympy

def find_fixed_point_coupling():
    """
    Calculates and explains the leading order expression for the Wilson-Fisher 
    fixed point coupling u* in phi^4 theory near d=4 dimensions.
    """
    # Define the symbols for the coupling 'u' and the small parameter 'epsilon'
    u, epsilon = sympy.symbols('u epsilon')

    print("This script finds the leading order expression for the fixed point coupling u*.")
    print("The theory is the phi^4 theory for a single scalar field near 4 dimensions (d = 4 - epsilon).")
    print("-" * 70)

    # Step 1: Define the one-loop beta function
    # The beta function describes how the coupling 'u' changes with the energy scale.
    # In a common normalization scheme, the one-loop beta function is:
    # β(u) = -εu + 3u²
    beta_function = -epsilon * u + 3 * u**2
    
    print("1. The one-loop Renormalization Group (RG) beta function for the coupling u is:")
    sympy.pprint(sympy.Eq(sympy.Symbol('β(u)'), beta_function), use_unicode=True)
    print("\n")

    # Step 2: Find fixed points by solving β(u*) = 0
    print("2. Fixed points (u*) are the values of u where the RG flow stops, i.e., where β(u*) = 0.")
    print("   We solve the equation:")
    equation = sympy.Eq(beta_function, 0)
    sympy.pprint(equation, use_unicode=True)
    print("\n")

    # Step 3: Solve the equation for u
    fixed_points = sympy.solve(equation, u)
    
    print(f"3. The solutions to this equation are u = {fixed_points[0]} and u = {fixed_points[1]}.")
    print(f"   - The solution u = {fixed_points[0]} is the 'Gaussian' or trivial fixed point, describing a non-interacting theory.")
    print(f"   - The non-zero solution is the non-trivial 'Wilson-Fisher' fixed point.\n")

    # Step 4: Display the final expression for the non-trivial fixed point
    wilson_fisher_fp = fixed_points[1]
    
    print("-" * 70)
    print("The leading order expression for the Wilson-Fisher fixed point coupling u* is:")
    
    # Final output formatted as an equation string
    final_equation = sympy.Eq(sympy.Symbol('u*'), wilson_fisher_fp)
    sympy.pprint(final_equation, use_unicode=True)
    
    # Fulfilling the request to output each number/symbol in the final equation
    print("\nIn this expression, the individual symbols and numbers are:")
    print(f"  - The fixed point coupling: u*")
    print(f"  - The small parameter: {epsilon}")
    print(f"  - A numerical constant: 3")


if __name__ == "__main__":
    find_fixed_point_coupling()