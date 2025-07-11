def solve_matrix_problem():
    """
    This function finds the number of natural numbers n for which there exist
    n real n-by-n matrices A_1, ..., A_n such that for any nonzero x in R^n,
    the vectors A_1x, ..., A_nx are linearly independent.

    The condition is equivalent to n <= rho(n), where rho(n) is the
    Radon-Hurwitz number. This function computes rho(n) for a range of n
    and identifies the values that satisfy the inequality.
    """

    def get_radon_hurwitz_number(n):
        """Calculates the Radon-Hurwitz number rho(n)."""
        if n <= 0:
            return 0
        
        # Decompose n = k * 2^m where k is odd.
        k = n
        m = 0
        while k % 2 == 0:
            k //= 2
            m += 1
            
        # Decompose m = 4d + c where 0 <= c <= 3.
        d = m // 4
        c = m % 4
        
        # Return rho(n) = 8d + 2^c
        return 8 * d + 2**c

    solutions = []
    # A theoretical analysis shows n > 8 cannot be a solution.
    # We will search up to n=32 to be thorough.
    for n in range(1, 33):
        if n <= get_radon_hurwitz_number(n):
            solutions.append(n)
            
    print("The natural numbers n for which the condition holds are:")
    for n_val in solutions:
        print(n_val)

    print("\nFor how many natural numbers n do there exist such matrices?")
    print(len(solutions))


if __name__ == "__main__":
    solve_matrix_problem()