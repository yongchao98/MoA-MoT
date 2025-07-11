def find_n_values():
    """
    This function finds the natural numbers n for which there exist n real n-by-n matrices
    A_1, ..., A_n such that for any nonzero vector x, the vectors A_1*x, ..., A_n*x are
    linearly independent.

    This condition is equivalent to n <= rho(n), where rho(n) is the Radon-Hurwitz number.
    """

    def get_rho(n):
        """Calculates the Radon-Hurwitz number rho(n)."""
        if n <= 0:
            return 0

        # Step 1: Write n = u * 2^m, where u is odd
        m = 0
        u = n
        while u % 2 == 0 and u != 0:
            u //= 2
            m += 1

        # Step 2: Write m = 4d + c
        d = m // 4
        c = m % 4

        # Step 3: Calculate rho(n) = 2^c + 8d
        rho_n = (2**c) + (8 * d)
        return rho_n

    solutions = []
    print("Checking for which natural numbers n the condition n <= rho(n) holds.\n")
    # We only need to check up to n=16 to be sure, as the inequality n <= rho(n)
    # diverges quickly for n > 8.
    for n in range(1, 17):
        rho_n = get_rho(n)
        if n <= rho_n:
            solutions.append(n)
            print(f"For n = {n}:")
            print(f"  - The Radon-Hurwitz number is rho({n}) = {rho_n}.")
            print(f"  - Checking the final equation: {n} <= {rho_n} is TRUE.")
            print(f"  - Therefore, n = {n} is a solution.\n")
        else:
            print(f"For n = {n}:")
            print(f"  - The Radon-Hurwitz number is rho({n}) = {rho_n}.")
            print(f"  - Checking the final equation: {n} <= {rho_n} is FALSE.\n")

    print("The possible values for n are:", solutions)
    print(f"In total, there are {len(solutions)} such natural numbers.")


find_n_values()