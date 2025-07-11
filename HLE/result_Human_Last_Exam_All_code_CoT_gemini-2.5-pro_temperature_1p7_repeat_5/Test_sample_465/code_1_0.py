def solve_chromatic_roots_problem():
    """
    This function analyzes the properties of chromatic and orbital chromatic roots
    and determines which of the given statements are true.
    """
    
    # Analysis of each statement
    # A: True. The set of real orbital chromatic roots is a subset of the set of real chromatic roots.
    # B: True. Chromatic roots can be complex numbers (e.g., for the Petersen graph).
    # C: True. Non-planar graphs can have negative real chromatic roots.
    # D: True. Non-integer real roots exist, like the golden ratio plus one.
    # E: False. It is a proven theorem that no chromatic roots exist in the interval (0, 1).

    true_statements = ['A', 'B', 'C', 'D']
    
    # Sort the letters alphabetically and join them into a single string.
    answer = "".join(sorted(true_statements))
    
    print("The true statements are A, B, C, and D.")
    print("The final answer as a sorted string is:")
    print(answer)

solve_chromatic_roots_problem()