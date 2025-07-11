import math

def solve():
    """
    This function calculates the query complexity for the two specified regimes.
    """

    # --- Regime 1: N = 2^sqrt(L) ---
    # In this regime, L = (log N)^2.
    # We compare two algorithms: Radix Sort and Comparison Sort.
    #
    # Radix Sort Complexity: O(N * L / log N)
    # Substituting L: O(N * (log N)^2 / log N) = O(N * log N)
    #
    # Comparison Sort Complexity: O(N * log N * log L)
    # Substituting L: O(N * log N * log((log N)^2)) = O(N * log N * 2 * log(log N))
    #
    # Radix sort is more efficient.
    # We convert O(N * log N) to the (a,b,c) format.
    # The target form is sqrt(N^a * (log N)^b * (log log N)^c).
    # We need to solve: N * log N = sqrt(N^a * (log N)^b * (log log N)^c)
    # Squaring both sides: N^2 * (log N)^2 = N^a * (log N)^b * (log log N)^c
    # By comparing exponents, we get:
    a1 = 2
    b1 = 2
    c1 = 0
    regime1_result = (a1, b1, c1)

    # --- Regime 2: N = 2^((log L)^2) ---
    # In this regime, L = 2^sqrt(log N).
    # We compare the same two algorithms.
    #
    # Radix Sort Complexity: O(N * L / log N)
    # Substituting L: O(N * 2^sqrt(log N) / log N). This grows super-polynomially in log N.
    #
    # Comparison Sort Complexity: O(N * log N * log L)
    # Substituting L: O(N * log N * log(2^sqrt(log N))) = O(N * log N * sqrt(log N))
    # This simplifies to O(N * (log N)^1.5).
    #
    # Comparison sort is more efficient.
    # We convert O(N * (log N)^1.5) to the (a,b,c) format.
    # We need to solve: N * (log N)^1.5 = sqrt(N^a * (log N)^b * (log log N)^c)
    # Squaring both sides: N^2 * (log N)^3 = N^a * (log N)^b * (log log N)^c
    # By comparing exponents, we get:
    a2 = 2
    b2 = 3
    c2 = 0
    regime2_result = (a2, b2, c2)

    # Print the final answer in the required format
    final_answer = f"{regime1_result},{regime2_result}"
    # The problem asks to output each number in the final equation.
    # For Regime 1: a=2, b=2, c=0
    # For Regime 2: a=2, b=3, c=0
    print(f"({a1},{b1},{c1}),({a2},{b2},{c2})")

solve()
<<<((2,2,0),(2,3,0))>>>