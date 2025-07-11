def hurwitz_radon(n):
    """
    Calculates the value of the Hurwitz-Radon function rho(n).

    The problem is equivalent to finding n for which there exists an n-dimensional
    subspace of n-by-n real matrices where every non-zero matrix is invertible.
    The maximum size of such a subspace is given by the Hurwitz-Radon function, rho(n).
    We need to find all natural numbers n such that rho(n) >= n.

    Formula for rho(n):
    1. Write n = u * 2^m, where u is odd.
    2. Write m = 4b + a, where a is in {0, 1, 2, 3}.
    3. rho(n) = 2^a + 8b.
    """
    if n <= 0:
        return 0
    
    # Step 1: Find u and m
    m = 0
    temp_n = n
    while temp_n > 0 and temp_n % 2 == 0:
        temp_n //= 2
        m += 1
    # u = temp_n, but we only need m

    # Step 2: Find a and b
    b = m // 4
    a = m % 4

    # Step 3: Calculate rho(n)
    rho_val = (2**a) + (8 * b)
    return rho_val

def solve_problem():
    """
    Finds the number of natural numbers n satisfying the condition.
    """
    print("Finding natural numbers n such that there exist n real n-by-n matrices A_1,...,A_n")
    print("such that for all nonzero x in R^n, A_1x,...,A_nx are linearly independent.")
    print("This is equivalent to finding n where the Hurwitz-Radon function rho(n) >= n.")
    print("-" * 30)

    solutions = []
    # We test n from 1 up to a reasonable limit to find the solutions.
    # The mathematical proof shows there are no solutions beyond n=8.
    limit = 32
    for n in range(1, limit + 1):
        rho_n = hurwitz_radon(n)
        is_solution = "Met" if rho_n >= n else "Not Met"
        if rho_n >= n:
            solutions.append(n)
            print(f"For n = {n:2d}: rho(n) = {rho_n:2d}. Condition rho(n) >= n is {is_solution}.  -> Solution found!")
        else:
            print(f"For n = {n:2d}: rho(n) = {rho_n:2d}. Condition rho(n) >= n is {is_solution}.")
    
    print("-" * 30)
    print(f"The possible values for n are: {solutions}")
    print(f"The number of such natural numbers is: {len(solutions)}")

if __name__ == '__main__':
    solve_problem()
