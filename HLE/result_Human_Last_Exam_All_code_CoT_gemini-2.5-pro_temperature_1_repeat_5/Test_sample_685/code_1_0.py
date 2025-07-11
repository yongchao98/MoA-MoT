def solve_nim_probability(n, m):
    """
    Calculates whether the first player in 2D-Generalized NIM has a winning
    position with a probability strictly more than 50% on a random n x m grid.

    This is equivalent to checking if the probability of a P-position is < 0.5.
    """
    print(f"Analyzing for n={n}, m={m}")

    # The condition is P(P-position) < 0.5, which is N_P / N_T < 1/2,
    # or 2 * N_P < N_T.

    # N_T is the total number of n x m binary matrices.
    N_T = 2**(n * m)

    # N_P is the number of P-positions.
    # A position is P if it's the zero matrix, or a rank-1 matrix uv^T where
    # u (n-dim) and v (m-dim) each have at least two 1s.

    # g(k) is the number of k-bit vectors with at least two 1s.
    g_n = 0
    if n >= 2:
        g_n = (2**n) - n - 1

    g_m = 0
    if m >= 2:
        g_m = (2**m) - m - 1

    # N_P = 1 (for the zero matrix) + number of rank-1 P-positions.
    N_P = 1 + g_n * g_m

    # The inequality to check is 2 * N_P < N_T.
    LHS = 2 * N_P
    RHS = N_T

    print("The final equation to determine the result is: 2 * N_P < N_T")
    print(f"We compare the left-hand side (LHS) and right-hand side (RHS):")
    # Output each number in the final equation.
    print(f"LHS = 2 * {N_P} = {LHS}")
    print(f"RHS = 2^({n}*{m}) = {RHS}")

    has_winning_chance = "strictly more than 50%" if LHS < RHS else "not strictly more than 50%"

    print(f"\nResult: The first player has a winning position with a probability {has_winning_chance}.")
    
    # The function f(n, m) would return 1 if LHS < RHS, and 0 otherwise.
    f_nm = 1 if LHS < RHS else 0
    print(f"f({n}, {m}) = {f_nm}")


# You can change these values to test different grid sizes.
n_val = 3
m_val = 4
solve_nim_probability(n_val, m_val)

print("\n" + "="*40 + "\n")

# Another example: the only case where the probability is not > 50%
n_val = 1
m_val = 1
solve_nim_probability(n_val, m_val)
