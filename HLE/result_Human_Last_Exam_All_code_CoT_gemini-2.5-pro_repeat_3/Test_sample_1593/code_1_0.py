def solve_complexity():
    """
    This function determines the query complexity for sorting bitstrings
    in two different regimes and prints the result in (a,b,c) notation.

    The complexity is determined by comparing two primary algorithms:
    1. Comparison Sort: O(N log N)
    2. Radix Sort with optimal block size: O(N * L / log N)

    The final complexity is Theta(sqrt(N^a * (log N)^b * (log log N)^c)).
    """

    # --- Regime 1: N = 2^sqrt(L) ---
    # In this regime, L = (log N)^2.
    # Cost of Comparison Sort: O(N log N)
    # Cost of Radix Sort: O(N * (log N)^2 / log N) = O(N log N)
    # The minimum complexity is O(N log N).
    # To find (a, b, c) for O(N log N):
    # N log N = sqrt(N^2 * (log N)^2) = sqrt(N^2 * (log N)^2 * (log log N)^0)
    # So, a=2, b=2, c=0.
    regime1_a = 2
    regime1_b = 2
    regime1_c = 0

    # --- Regime 2: N = 2^((log L)^2) ---
    # In this regime, L = 2^sqrt(log N).
    # Cost of Comparison Sort: O(N log N)
    # Cost of Radix Sort: O(N * 2^sqrt(log N) / log N)
    # Since 2^sqrt(log N) grows much faster than (log N)^2, the Radix Sort
    # complexity is higher than the Comparison Sort complexity.
    # The minimum complexity is O(N log N).
    # The (a,b,c) representation is the same as in Regime 1.
    regime2_a = 2
    regime2_b = 2
    regime2_c = 0

    # The problem asks for the answer in the format (a1,b1,c1),(a2,b2,c2)
    # where each number is part of the final equation for complexity.
    final_answer = f"({regime1_a},{regime1_b},{regime1_c}),({regime2_a},{regime2_b},{regime2_c})"
    print(final_answer)

solve_complexity()