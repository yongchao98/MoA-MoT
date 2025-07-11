def solve():
    """
    Calculates the sum based on the derived formula and prints the steps.
    """
    # The size of the vocabulary, as given in the problem.
    n = 99

    # The mathematical derivation shows that the sum S is equal to (n+1)^(n-1).
    # We will calculate this value step by step.

    print("The problem is to calculate the sum S = sum(a(w)) for all w in V^n.")
    print(f"The size of the vocabulary is n = {n}.")
    print("Through mathematical derivation, the sum S simplifies to the formula (n+1)^(n-1).")
    print("\nNow, we substitute n = 99 into the formula to find the value:")

    # Step 1: Calculate the base of the power.
    base = n + 1
    # Step 2: Calculate the exponent of the power.
    exponent = n - 1

    print(f"S = ({n} + 1)^({n} - 1)")
    print(f"S = {base}^{exponent}")

    # Step 3: Express the result as a power of 10.
    # Since the base is 100, which is 10^2, we can rewrite the expression.
    final_exponent = 2 * exponent
    
    print("\nTo express this as a power of 10:")
    print(f"S = (10^2)^{exponent}")
    print(f"S = 10^(2 * {exponent})")
    print(f"S = 10^{final_exponent}")

solve()