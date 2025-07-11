import numpy as np

def illustrate_proof():
    """
    This function illustrates the theoretical proof with a concrete example.
    """
    # 1. Define an example Markov Chain that fits the criteria.
    # Consider a biased random walk on the non-negative integers {0, 1, 2, ...}.
    # Let the probability of moving right (p) be greater than moving left (q).
    # This chain is known to be transient (and thus not positive recurrent).
    p = 0.6
    q = 1.0 - p

    # The state space is Sigma = {0, 1, 2, ...}
    # Transitions: p(i, i+1) = p, p(i, i-1) = q for i >= 1.
    # At the boundary, let's say it's a reflecting barrier: p(0, 1) = 1.

    # 2. Define the finite set A and the function f.
    # Let A be the set containing only the origin.
    A = {0}
    # Let f(x) = (p/q)^x. Since p > q, (p/q) > 1, so f(x) -> infinity as x -> infinity.
    def f(x, p_val, q_val):
        return (p_val / q_val)**x

    # 3. Verify the submartingale condition for x not in A.
    # Condition: E[f(X_1) | X_0=x] - f(x) >= 0 for x > 0.
    # Drift at x: p*f(x+1) + q*f(x-1) - f(x)
    def calculate_drift(x, p_val, q_val):
        expected_f_next = p_val * f(x + 1, p_val, q_val) + q_val * f(x - 1, p_val, q_val)
        return expected_f_next - f(x, p_val, q_val)

    # 4. Walk through the logic of the contradiction.
    print("--- Illustrating the Proof with an Example ---")
    print(f"Consider a biased random walk with p={p}, q={q}.")
    print("The chain is on {0, 1, 2, ...} with a reflecting barrier at 0.")
    print(f"The finite set is A = {A}.")
    print(f"The function is f(x) = (p/q)^x = ({p/q:.2f})^x.")
    print("This function satisfies f(x) -> infinity as x -> infinity.")
    
    # Check drift condition
    x_check = 5
    drift = calculate_drift(x_check, p, q)
    print(f"The drift condition E[f(X_1)|X_0=x] - f(x) >= 0 for x not in A is satisfied.")
    print(f"For x={x_check}, the drift is {drift:.4f}, which is positive.\n")

    print("The proof by contradiction proceeds as follows:")
    print("1. Assume the chain is positive recurrent (and therefore recurrent).")
    print("2. This implies that for any starting state x, f(x) <= max_{y in A} f(y).")

    # The "final equation" is the contradiction itself.
    max_f_on_A = f(list(A)[0], p, q)
    print(f"3. In our example, max_{{y in A}} f(y) = f(0) = {max_f_on_A:.4f}.")
    print(f"4. So, the assumption implies that f(x) <= {max_f_on_A:.4f} for ALL states x.")

    # 5. Show the contradiction.
    x_test = 2
    f_x_test = f(x_test, p, q)
    print("\nLet's test this implication for state x = 2.")
    print("The derived inequality from the recurrence assumption is:")
    print(f"f({x_test}) <= {max_f_on_A}")
    print("The actual values are:")
    # This is the "equation" from the prompt, showing the contradiction
    print(f"{f_x_test:.4f} <= {max_f_on_A}")
    print(f"This is FALSE.")
    
    print("\n--- Conclusion ---")
    print("The contradiction shows the initial assumption was wrong.")
    print("Therefore, the Markov chain cannot be positive recurrent.")

if __name__ == '__main__':
    illustrate_proof()
    print("\nFinal Answer: Yes, one can conclude the chain is not positive recurrent.")