def solve_matrix_problem():
    """
    This function solves the problem by leveraging a known theorem in mathematics.

    The problem is equivalent to finding the number of natural numbers 'n' for which
    the sphere S^(n-1) is parallelizable. A theorem by Adams, Bott, Milnor, and
    Kervaire states that this occurs only for n = 1, 2, 4, and 8.
    """

    # The set of natural numbers n for which the condition holds.
    solution_set = [1, 2, 4, 8]

    # The problem asks for "how many" such natural numbers exist.
    number_of_solutions = len(solution_set)

    # To satisfy the prompt "output each number in the final equation",
    # we first display the set of solutions and then its size.
    
    # Create a string representation of the set, e.g., "{1, 2, 4, 8}"
    set_string = "{" + ", ".join(map(str, solution_set)) + "}"
    
    print(f"The set of natural numbers n for which the property holds is S = {set_string}.")
    print(f"The number of such natural numbers is the size of the set S, which is |S| = {number_of_solutions}.")

solve_matrix_problem()