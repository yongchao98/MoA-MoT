import math

def solve():
    """
    Solves for the query complexity in the (a,b,c) notation for the two regimes.
    """

    # The (a,b,c) notation represents Theta(sqrt(N^a * (log N)^b * (log log N)^c)).
    # We need to find the complexity Q and then find a, b, c such that
    # Q^2 is proportional to N^a * (log N)^b * (log log N)^c.

    results = []

    # --- Regime 1: N = 2^sqrt(L) ---
    # This means log2(N) = sqrt(L), so L = (log2 N)^2.

    # We compare the complexity of two algorithms:
    # 1. Comparison Sort: Q1 = O(N log N)
    # 2. Hash-and-Sort: Q2 = O(N + U log U) where U = min(N, 2^L)

    # To find U, we compare N and 2^L by comparing log2(N) and L.
    # log2(N) vs L is log2(N) vs (log2 N)^2.
    # For N > 4 (i.e., log2(N) > 2), we have (log2 N)^2 > log2 N, so L > log2 N.
    # This implies 2^L > N.
    # Therefore, U = min(N, 2^L) = N.

    # The complexity of Hash-and-Sort is O(N + N log N) = O(N log N).
    # The minimum of the two algorithms is Theta(N log N).

    # Convert Theta(N log N) to (a,b,c) notation.
    # Q = N * log(N)
    # Q^2 = (N * log(N))^2 = N^2 * (log N)^2
    # Q^2 = N^2 * (log N)^2 * (log log N)^0
    # So, a=2, b=2, c=0.
    regime1_abc = (2, 2, 0)
    results.append(regime1_abc)

    # --- Regime 2: N = 2^((log2 L)^2) ---
    # This means log2 N = (log2 L)^2.

    # Again, we compare the two algorithms. We need to find U = min(N, 2^L).
    # We compare log2(N) and L.
    # log2(N) vs L is (log2 L)^2 vs L.
    # Asymptotically, for large L, L grows faster than (log2 L)^2.
    # So, we assume L > log2 N for the general case in this regime.
    # This implies 2^L > N.
    # Therefore, U = min(N, 2^L) = N.

    # The complexity of Hash-and-Sort is O(N + N log N) = O(N log N).
    # The minimum of the two algorithms is again Theta(N log N).

    # Convert Theta(N log N) to (a,b,c) notation.
    # This is the same as in Regime 1.
    # Q = N * log(N)
    # Q^2 = N^2 * (log N)^2
    # a=2, b=2, c=0.
    regime2_abc = (2, 2, 0)
    results.append(regime2_abc)

    # Format the final answer string
    final_answer = ",".join([f"({a},{b},{c})" for a, b, c in results])
    print(final_answer)

solve()