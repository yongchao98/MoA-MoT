def find_valid_n_values():
    """
    Finds all natural numbers n that satisfy the condition from the problem.

    The condition is equivalent to n <= rho(n), where rho(n) is the
    Radon-Hurwitz number.
    """

    def get_rho(n):
        """Calculates the Radon-Hurwitz number rho(n)."""
        if n <= 0:
            return 0
        
        # Write n as 2^a * b where b is odd
        a = 0
        b = n
        while b % 2 == 0:
            b //= 2
            a += 1
        
        # Write a as 4q + r
        q = a // 4
        r = a % 4
        
        # rho(n) = 2^r + 8q
        return 2**r + 8*q

    solutions = []
    # Based on the analysis of the inequality n <= rho(n), solutions
    # are expected to be small. We check up to n=20, which is sufficient
    # to find all solutions and show that they stop appearing.
    print("Checking condition n <= rho(n) for n = 1 to 20:")
    for n in range(1, 21):
        rho_n = get_rho(n)
        if n <= rho_n:
            status = "holds"
            solutions.append(n)
        else:
            status = "does not hold"
        print(f"For n = {n:2d}, rho({n}) = {rho_n:2d}. Condition {n} <= {rho_n} {status}.")

    print("\nThe natural numbers n for which the condition holds are:", solutions)
    print("The total count of such natural numbers is:", len(solutions))

# Run the function to find and print the results.
find_valid_n_values()