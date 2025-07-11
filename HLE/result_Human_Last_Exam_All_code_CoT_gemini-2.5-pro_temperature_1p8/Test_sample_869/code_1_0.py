import math

def solve_probability():
    """
    Solves the given combinatorial probability problem.
    """
    # Problem parameters
    N = 25  # Total items
    T = 5   # Number of types / individuals
    K = 5   # Items per type / per individual

    # Step 1: Calculate the total number of possible distributions (S)
    # This is the number of ways to arrange N items with T groups of K identical items.
    # S = N! / (K!)^T
    S = math.factorial(N) // (math.factorial(K)**T)

    # Step 2: Calculate the number of favorable distributions (F)
    # This involves finding valid distribution matrices C, counting how many ways
    # each can be formed, and summing them up.

    # We found three types of matrices for a fixed specialty mapping (e.g., person i -> type i).

    # Case 1: C = 5*I. Each person gets 5 items of their specialty type.
    # Number of such matrices for a fixed mapping: 1
    # Arrangements for this matrix: (K!)^T / (K!)^T = 1
    num_matrices_c1 = 1
    ways_c1 = 1
    f_case1 = num_matrices_c1 * ways_c1

    # Case 2: C = 4*I + P, where P is a derangement permutation matrix.
    # Number of derangements of 5 items: D(5) = 44
    num_matrices_c2 = 44
    # Arrangements for one such matrix: (K!)^T / ((K-1)!^T * 1!^T) = K^T
    ways_c2 = K**T
    f_case2 = num_matrices_c2 * ways_c2

    # Case 3: C = 3*I + P + P^T, where P is a 5-cycle permutation matrix.
    # Number of 5-cycles on 5 vertices: (T-1)! / 2 = 12
    num_matrices_c3 = math.factorial(T - 1) // 2
    # Arrangements for one such matrix: (K!)^T / ((K-2)!^T * 1!^(2T)) = (K*(K-1))^T
    ways_c3 = (K * (K - 1))**T
    f_case3 = num_matrices_c3 * ways_c3

    # The number of favorable arrangements for a *fixed* specialty mapping
    F_fixed_mapping = f_case1 + f_case2 + f_case3

    # Total number of favorable arrangements F is F_fixed_mapping multiplied by the
    # number of possible specialty mappings (T!).
    F = math.factorial(T) * F_fixed_mapping

    # Step 3: Calculate the probability P = F / S
    P = F / S

    # Final Output
    print("This problem calculates the probability of a specific type of distribution of items.")
    print("The calculation proceeds by finding the total number of distributions (S) and")
    print("the number of 'favorable' distributions (F) that meet the problem's criteria.")
    print("\n--- Calculation Details ---")
    print(f"The total number of ways to distribute the items is S:")
    print(f"S = {N}! / ({K}!)**{T}")
    print(f"S = {S}")

    print(f"\nThe number of favorable distributions is F:")
    print(f"This is based on summing the ways for three types of valid distribution patterns,")
    print(f"multiplied by the {T}! ways to assign specialties to individuals.")
    f_calc_str = (
        f"F = {T}! * ({num_matrices_c1}*{ways_c1} + {num_matrices_c2}*{K}**{T} + "
        f"{num_matrices_c3}*({K}*({K}-1))**{T})"
    )
    print(f_calc_str)
    f_sum_str = (
        f"F = {math.factorial(T)} * ({f_case1} + {f_case2} + {f_case3})"
    )
    print(f_sum_str)
    print(f"F = {math.factorial(T)} * ({F_fixed_mapping})")
    print(f"F = {F}")


    print("\n--- Final Probability Equation ---")
    print("The probability P is the ratio F / S.")
    print(f"P = {F} / {S}")
    print(f"P ≈ {P:.12f}")
    
    return P

# Execute the solution
final_probability = solve_probability()
# The final answer is wrapped according to the format.
# <<<final_probability>>> will not be displayed in the output.