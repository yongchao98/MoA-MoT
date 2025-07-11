import math

def solve():
    """
    Solves the query complexity problem for two regimes.
    """

    print("### Step 1: Establish General Complexity Formulas ###")
    print("We have two primary algorithmic strategies to sort N bitstrings of length L:")
    print("1. Comparison-based Sort: Treat each string as an atomic item and sort. This can be done by first finding unique strings with N 'H' queries and then sorting the M unique strings with O(M log M) 'C' queries. In the worst case (M=N), this costs O(N log N) queries.")
    print("2. Radix Sort: Process strings in chunks. An optimized radix sort algorithm has a query complexity of O(L*N / log N).")
    print("\nThe overall query complexity Q is the minimum of these two approaches:")
    print("Q = min(O(N log N), O(L*N / log N))")
    print("-" * 20)

    # --- Regime 1 ---
    print("\n### Step 2: Analyze Regime 1: N = 2^sqrt(L) ###")
    print("First, we express L in terms of N:")
    print("log2(N) = sqrt(L)  =>  L = (log2(N))^2")
    print("\nNow, we evaluate the complexity of the Radix Sort approach in this regime:")
    print("O(L*N / log N) = O((log N)^2 * N / log N) = O(N log N)")
    print("\nIn this regime, both algorithms yield a complexity of O(N log N).")
    print("Therefore, the query complexity is Q1 = Theta(N log N).")
    print("-" * 20)

    # --- Regime 2 ---
    print("\n### Step 3: Analyze Regime 2: N = 2^((log2(L))^2) ###")
    print("The relationship is log2(N) = (log2(L))^2.")
    print("We need to compare the two complexities, which means comparing 'log N' and 'L / log N'.")
    print("This is equivalent to comparing (log N)^2 and L.")
    print("(log N)^2 = ((log2(L))^2)^2 = (log2(L))^4.")
    print("For large L, L grows faster than any power of log(L). So, L > (log2(L))^4.")
    print("This means L > (log N)^2, which implies L / log N > log N.")
    print("Therefore, the O(N log N) complexity from the comparison-based sort is lower.")
    print("The query complexity is Q2 = Theta(N log N).")
    print("-" * 20)

    # --- (a,b,c) notation ---
    print("\n### Step 4: Convert to (a,b,c) Notation ###")
    print("Both regimes have a query complexity of Q = Theta(N log N).")
    print("We need to find a, b, c such that sqrt(N^a * (log N)^b * (log log N)^c) = Q.")
    print("Squaring both sides of the equation:")
    print("N^a * (log N)^b * (log log N)^c = Q^2 = (N log N)^2 = N^2 * (log N)^2")
    print("\nBy comparing the exponents on both sides, we get the final equation values:")
    a = 2
    b = 2
    c = 0
    print(f"a = {a}")
    print(f"b = {b}")
    print(f"c = {c}")
    print("\nThis gives the complexity class (2, 2, 0) for both regimes.")
    print("-" * 20)

    # --- Final Answer ---
    print("\n### Step 5: Final Answer ###")
    result_regime1 = (a, b, c)
    result_regime2 = (a, b, c)
    final_answer = f"({result_regime1[0]},{result_regime1[1]},{result_regime1[2]}),({result_regime2[0]},{result_regime2[1]},{result_regime2[2]})"
    print("The final answer in the format (a,b,c),(a,b,c) is:")
    print(final_answer)

if __name__ == '__main__':
    solve()