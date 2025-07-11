def solve_markov_chain_property():
    """
    Illustrates the property of a non-positive recurrent Markov chain.

    We consider a biased random walk on the non-negative integers Z+ = {0, 1, 2, ...}.
    - State Space Sigma: Z+
    - Transition probabilities:
        - For x > 0: p(x, x+1) = p, p(x, x-1) = q, with p > 1/2.
        - For x = 0: p(0, 1) = 1 (reflecting barrier).
    This chain is irreducible. Since p > 1/2, it has a drift to +infinity and is transient.

    - Finite set A = {0}.
    - Function f(x) = (p/q)^x. Since p > q, (p/q) > 1 and f(x) -> infinity as x -> infinity.

    We will verify the condition: E[f(X_1) | X_0=x] - f(x) >= 0 for x not in A (i.e., x > 0).
    """
    # 1. Define the parameters for our example chain.
    p = 0.6
    q = 1 - p

    # 2. Define the function f(x)
    ratio = p / q
    def f(x):
        return ratio**x

    # 3. Choose a state x outside of A to test the condition.
    # A = {0}, so let's choose x = 2.
    x_test = 2
    
    # We need to verify that p*f(x+1) + q*f(x-1) - f(x) >= 0
    # For x_test = 2, the next states are x+1=3 and x-1=1.
    fx = f(x_test)
    fx_plus_1 = f(x_test + 1)
    fx_minus_1 = f(x_test - 1)

    # 4. Calculate the drift Lf(x) = E[f(X_1) | X_0=x] - f(x)
    drift = p * fx_plus_1 + q * fx_minus_1 - fx

    # 5. Print the verification and conclusion.
    print("Let's verify the condition for our example Markov Chain.")
    print(f"Parameters: p = {p}, q = {q}. The set A = {{0}}.")
    print(f"The function is f(x) = ({p}/{q})^x = {ratio:.2f}^x.")
    print(f"\nWe test the condition at state x = {x_test} (which is not in A).")
    print("The condition is: p*f(x+1) + q*f(x-1) - f(x) >= 0")
    
    # Output the equation with numbers, as requested.
    print("\nThe final equation with numbers is:")
    print(f"{p} * f({x_test+1}) + {q} * f({x_test-1}) - f({x_test}) = ", end="")
    print(f"{p} * {fx_plus_1:.4f} + {q} * {fx_minus_1:.4f} - {fx:.4f} = {drift:.4f}")

    print(f"\nThe result {drift:.4f} is indeed non-negative, so the condition holds.")
    print("\nAs shown by the theoretical argument and illustrated here, a chain with these properties cannot be positive recurrent.")

solve_markov_chain_property()