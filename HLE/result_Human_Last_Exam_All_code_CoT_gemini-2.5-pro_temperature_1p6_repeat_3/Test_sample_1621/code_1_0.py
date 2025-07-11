def solve_matrix_problem():
    """
    This function determines for how many natural numbers n there exist n real n-by-n matrices
    A_1, ..., A_n such that for all nonzero x in R^n, the vectors A_1*x, ..., A_n*x are
    linearly independent.

    The problem can be solved by analyzing the properties of the determinant of the matrix
    formed by the resulting vectors. This leads to a connection with deep results in
    algebraic topology, specifically the parallelizability of spheres.

    The solution is that such matrices exist only for n = 1, 2, 4, and 8. These values
    correspond to the dimensions of the real normed division algebras (reals, complex
    numbers, quaternions, and octonions).

    Therefore, there are 4 such natural numbers.
    """

    # The values of n for which the condition holds, based on the mathematical proof.
    valid_n_values = [1, 2, 4, 8]

    # The count of such numbers.
    count = len(valid_n_values)

    print("The natural numbers n for which the condition holds are:")
    # We are asked to output each number in the final equation.
    # The 'equation' here is the set of numbers that satisfy the condition.
    for number in valid_n_values:
        print(number)

    print("\nThe total number of such natural numbers is:")
    print(count)

solve_matrix_problem()