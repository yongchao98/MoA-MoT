import sympy

def find_fixed_point_coupling():
    """
    Calculates and explains the leading order expression for the fixed point
    coupling u* in phi^4 theory near d=4 dimensions (d=4-epsilon).
    """

    # Define the symbols for our equation
    u_star = sympy.Symbol('u*')
    epsilon = sympy.Symbol('ε')
    pi = sympy.pi

    print("Step 1: The Beta Function for ϕ⁴ Theory")
    print("------------------------------------------")
    print("In d = 4 - ε dimensions, the one-loop beta function for the coupling u is:")
    print("β(u) = -ε*u + (3 / (16*π²)) * u²\n")

    # Define the beta function expression using sympy
    # The coefficients are standard results from quantum field theory.
    beta_function = -epsilon * u_star + (3 / (16 * pi**2)) * u_star**2

    print("Step 2: The Fixed Point Condition")
    print("-----------------------------------")
    print("A fixed point u* is found where the beta function is zero:")
    print("β(u*) = 0")
    print(f"{sympy.pretty(beta_function)} = 0\n")

    print("Step 3: Solving for the Fixed Point u*")
    print("---------------------------------------")
    # Solve the equation beta(u*) = 0 for u*
    solutions = sympy.solve(beta_function, u_star)

    # We are interested in the non-trivial (Wilson-Fisher) fixed point
    wilson_fisher_fp = None
    for sol in solutions:
        if sol != 0:
            wilson_fisher_fp = sol
            break

    print("The equation has two solutions:")
    print(f"  1. Trivial (Gaussian) fixed point: u* = {solutions[0]}")
    print(f"  2. Non-trivial (Wilson-Fisher) fixed point: u* = {solutions[1]}\n")

    print("Step 4: The Leading Order Expression")
    print("------------------------------------")
    print("The leading order expression for the non-trivial fixed point is derived from:")
    print("-ε + (3 / (16*π²)) * u* = 0")
    print("Which gives:")
    
    # Extracting the numbers for the final formatted output
    num = 16
    den = 3
    
    print(f"u* = ({num} * π² / {den}) * ε")


# Execute the function to print the solution steps
find_fixed_point_coupling()