def solve_max_queens():
    """
    Calculates the maximum number m such that m white queens and m black queens
    can coexist on an N x N chessboard without same-color attacks.
    This is demonstrated for N=16.
    """
    N = 16

    # This problem is equivalent to finding if two disjoint N-queens solutions exist.
    # If they do, m_max = N. We will construct and verify a pair of solutions.

    # A queen configuration can be represented by a permutation array `p` of length N,
    # where p[i] is the column of the queen in row `i`. We use 0-based indexing.

    # Step 1: Construct the first solution (P1) using a known formula for N where N mod 6 = 4.
    # The permutation is p[i] = 2i+1 for i in [0, N/2-1] and p[i] = 2(i-N/2) for i in [N/2, N-1].
    p1 = [0] * N
    half_n = N // 2
    for i in range(half_n):
        p1[i] = 2 * i + 1
    for i in range(half_n, N):
        p1[i] = 2 * (i - half_n)

    # Step 2: Construct the second solution (P2) by creating a horizontal reflection of P1.
    # A reflection of a valid solution is also a valid solution.
    p2 = p1[::-1]

    def is_valid_solution(p, n):
        """Checks if a permutation represents a valid N-queens solution."""
        for i in range(n):
            for j in range(i + 1, n):
                # Check for diagonal attacks. Row/column attacks are impossible by representation.
                if abs(p[i] - p[j]) == abs(i - j):
                    return False
        return True

    # Step 3: Verify both solutions.
    p1_is_valid = is_valid_solution(p1, N)
    p2_is_valid = is_valid_solution(p2, N)

    # Step 4: Verify that the two solutions are disjoint.
    # This means that for any given row `i`, the queens are in different columns.
    # i.e., p1[i] != p2[i] for all i.
    are_disjoint = True
    for i in range(N):
        if p1[i] == p2[i]:
            are_disjoint = False
            break

    # Step 5: Conclude and print the result.
    # If we found two valid, disjoint N-queens solutions, then m=N is possible.
    # Since m cannot be > N, this is the maximum value.
    if p1_is_valid and p2_is_valid and are_disjoint:
        max_m = N
        print(f"A valid 16-queens solution (P1) is: {[x + 1 for x in p1]}")
        print(f"A second valid 16-queens solution (P2) is: {[x + 1 for x in p2]}")
        print("The two solutions are disjoint, proving that 16 white and 16 black queens can coexist.")
        print(f"The maximum number m is: {max_m}")
    else:
        # This fallback should not be reached with the chosen construction.
        # It would imply the construction is flawed.
        print("Could not constructively verify the solution. Based on the simpler interpretation (2m <= 16), m would be 8.")

solve_max_queens()