import sympy

def solve_graph_problem():
    """
    Analyzes the system of equations derived from the problem's constraints
    to determine if a solution is possible.
    """
    # Define the variables as symbolic integers
    n, n0, n1, n2 = sympy.symbols('n n0 n1 n2', integer=True)

    # Equation 1: The sum of vertices in each category equals the total number of vertices, n.
    # n0 + n1 + n2 = n
    eq1 = sympy.Eq(n0 + n1 + n2, n)
    print("Equation 1 (Total vertices):")
    print(f"{eq1.lhs} = {eq1.rhs}")
    print("-" * 30)

    # Equation 2: The double-counting argument for vertex-C5 incidences.
    # Sum of c(v) = (0*n0) + (1*n1) + (2*n2) = 5n
    eq2 = sympy.Eq(n1 + 2 * n2, 5 * n)
    print("Equation 2 (Double-counting C5s):")
    print(f"{eq2.lhs} = {eq2.rhs}")
    print("-" * 30)
    
    # We solve this system of 2 equations for the 3 variables (n0, n1, n2).
    # We can express n0 and n1 in terms of n and n2.
    solution = sympy.solve([eq1, eq2], [n0, n1])
    
    print("Solving for n0 and n1 in terms of n and n2 gives:")
    n0_expr = solution[n0]
    n1_expr = solution[n1]
    
    print(f"n0 = {n0_expr}")
    print(f"n1 = {n1_expr}")
    print("-" * 30)
    
    # Now, we apply the physical constraints that these counts cannot be negative.
    # n0, n1, n2 must be >= 0.
    # Also, n2 (number of vertices belonging to 2 C5s) cannot exceed the total number of vertices n.
    
    print("Applying physical constraints:")
    print("Constraint 1: n0 must be non-negative (n0 >= 0).")
    # From n0 = n2 - 4*n, this gives:
    # n2 - 4*n >= 0  =>  n2 >= 4*n
    print(f"From n0 = {n0_expr}, we get the inequality: n2 >= 4*n")

    print("\nConstraint 2: n2 must be less than or equal to the total number of vertices n (n2 <= n).")
    
    print("\nCombining these constraints, we get a final inequality for n:")
    # We have 4*n <= n2 and n2 <= n.
    # This implies 4*n <= n.
    print(f"4*n <= n2 <= n  =>  4*n <= n")
    print("The final inequality is 4*n <= n, which simplifies to 3*n <= 0.")
    print("-" * 30)
    
    print("Conclusion:")
    print("The problem asks for a composite number n, which must be a positive integer (n > 0).")
    print("However, the inequality 3*n <= 0 can only be true for n <= 0.")
    print("This is a logical contradiction.")
    print("\nTherefore, no graph satisfying all the given properties exists for any positive n.")
    print("As a result, there is no smallest composite n for which such a graph exists.")

solve_graph_problem()