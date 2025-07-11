def solve_matrix_problem():
    """
    This function solves the problem of finding for how many natural numbers n
    there exist n real n-by-n matrices A_1, ..., A_n such that for all
    nonzero x in R^n, the vectors A_1x, ..., A_nx are linearly independent.
    
    The solution relies on a deep mathematical theorem by Frank Adams, which
    states that this is only possible for n in {1, 2, 4, 8}.
    Our explanation details why odd n > 1 are ruled out by a simpler argument.
    """
    
    # Based on the mathematical reasoning, the possible values for n are 1, 2, 4, and 8.
    # These correspond to the dimensions of the real division algebras.
    solution_n = [1, 2, 4, 8]
    
    print("The problem asks for how many natural numbers n do there exist n real n-by-n matrices")
    print("A_1, ..., A_n such that for all nonzero x in R^n, the vectors A_1x, ..., A_nx")
    print("are linearly independent.\n")
    
    print("The possible values for n are a result of a theorem by Frank Adams and are related")
    print("to the dimensions of real division algebras (Real Numbers, Complex Numbers, Quaternions, Octonions).")

    # The phrasing "output each number in the final equation" is interpreted
    # as listing the solutions and the final count.
    # For example, "The set of solutions n is {1, 2, 4, 8}, and the count is 4."
    
    print("\nThe set of solutions for n is:")
    print(f"n = {solution_n}")
    
    count = len(solution_n)
    print(f"\nThe number of such natural numbers n is: {count}")

solve_matrix_problem()