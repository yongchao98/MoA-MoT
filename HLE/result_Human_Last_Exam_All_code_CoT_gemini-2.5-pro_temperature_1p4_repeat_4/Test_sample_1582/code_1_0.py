def demonstrate_proof_by_contradiction():
    """
    This script explains the proof that the described Markov chain cannot be
    positive recurrent. It uses a set of hypothetical values to make the
    logical steps and the final contradiction clear.
    """

    # --- Setup of the hypothetical scenario based on the problem's premises ---

    # 1. Let A be the finite set from the problem. Since A is finite, f(x) has a
    #    maximum value on A. Let's assume a hypothetical value for this maximum.
    max_f_on_A = 50.0

    # 2. We choose a constant C0 that is greater than the maximum of f on A.
    C0 = 100.0

    # 3. We define a set B = {x in Sigma | f(x) <= C0}. Since f(x) -> infinity,
    #    this set B must be finite.

    # 4. We can always find a state x0 where f(x0) is larger than C0.
    #    Let's choose a starting state x0 with a specific value for f(x0).
    f_x0 = 150.0

    # --- The Argument ---

    print("Yes, one can conclude the chain is not positive recurrent.")
    print("The reasoning is a proof by contradiction, illustrated below with example values:")
    print("-" * 75)

    print(f"Step 1: Assume the chain IS positive recurrent (and therefore recurrent).")
    print(f"Step 2: Let's pick a constant C0 = {C0}, which is greater than the max of f(x) on set A (e.g., {max_f_on_A}).")
    print(f"Step 3: Let's pick a starting state x0 where f(x0) > C0. We'll use f(x0) = {f_x0}.")
    print(f"         Our choice means: f(x0) > C0")

    print("\nStep 4: Now, let's see what our assumption implies:")
    print(f"  - If the chain is recurrent, it must eventually hit the finite set B = {{x | f(x) <= {C0}}}.")
    print("  - The Optional Stopping Theorem applied to the submartingale-like process tells us:")
    print("    f(x0) <= E[f(X_tau_B)], where X_tau_B is the state where the chain first hits B.")
    print(f"  - By definition of B, f(X_tau_B) is always <= {C0}. So its expectation is also <= {C0}.")
    print(f"  - This leads to the conclusion: f(x0) <= {C0}")

    print("\nStep 5: The Contradiction.")
    print("Our logic has led to two conflicting statements which form our final equation:")

    # The final equation demonstrating the contradiction, with each number printed.
    print(f"    - From our choice in Step 3: {f_x0} > {C0}")
    print(f"    - From our proof in Step 4:  {f_x0} <= {C0}")

    print("\nThis is a logical impossibility. The initial assumption must be false.")
    print("Conclusion: The Markov chain is not positive recurrent.")

demonstrate_proof_by_contradiction()