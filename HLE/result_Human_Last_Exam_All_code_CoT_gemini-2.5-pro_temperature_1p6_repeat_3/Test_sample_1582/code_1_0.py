# This script outlines the logical proof by contradiction.

# --- Given Conditions ---
# 1. An irreducible Markov chain (let's call it X_n) on a state space Sigma.
# 2. A finite subset A of Sigma.
# 3. A non-negative function f: Sigma -> R_+ such that f(x) -> infinity as x -> infinity.
#    This implies that for any constant C, the set {x | f(x) <= C} is finite.
# 4. For all x not in A, E[f(X_1) | X_0 = x] >= f(x).

def prove_not_positive_recurrent():
    """
    This function lays out the proof that the Markov chain cannot be positive recurrent.
    """
    # Step 1: Assume the chain is positive recurrent. This is our hypothesis for contradiction.
    # A consequence is that the chain must also be recurrent.
    # For a recurrent irreducible chain, the first hitting time of any finite set is finite with probability 1.
    # Let tau_A be the first time the chain enters the set A.
    # So, P(tau_A < infinity) = 1.

    # Step 2: Use the property f(x) -> infinity.
    # Let M_A = max(f(y) for y in A). Since A is finite, M_A is a finite number.
    # Since f(x) -> infinity, there must exist states outside of A.
    # We can choose a starting state, x0, such that f(x0) > M_A.
    # By definition of M_A, this x0 cannot be in A.
    x0 = "a chosen state not in A"
    M_A = "max value of f on set A"
    print(f"Assertion 1: We select a state x0 such that f(x0) > {M_A}")

    # Step 3: Define a submartingale and apply the Optional Stopping Theorem (OST).
    # Let's define a new process Z_n = f(X_{n_and_tau_A}), where n_and_tau_A is min(n, tau_A).
    # The given conditions ensure Z_n is a non-negative submartingale.
    # The OST can be applied to Z_n and the stopping time tau_A.
    # The theorem states: E[Z_{tau_A} | X_0 = x0] >= Z_0
    # Substituting the definitions:
    # E[f(X_tau_A) | X_0 = x0] >= f(x0)
    LHS = "E[f(X_tau_A) | X_0 = x0]" # Left Hand Side of the inequality
    RHS = "f(x0)" # Right Hand Side
    print(f"Assertion 2 (from OST): {LHS} >= {RHS}")

    # Step 4: Analyze the left-hand side of the inequality from Step 3.
    # By definition of tau_A, the state X_tau_A must be in the set A.
    # Therefore, the value of the function f at this state is at most M_A.
    # f(X_tau_A) <= M_A
    # Taking the expectation of both sides:
    # E[f(X_tau_A) | X_0 = x0] <= M_A
    print(f"Assertion 3 (from properties of A): {LHS} <= {M_A}")

    # Step 5: Combine the results to reveal the contradiction.
    # From Assertion 2: f(x0) <= E[f(X_tau_A) | X_0 = x0]
    # From Assertion 3: E[f(X_tau_A) | X_0 = x0] <= M_A
    # Combining these gives: f(x0) <= M_A
    print(f"Conclusion from Assertions 2 & 3: f(x0) <= {M_A}")

    # The contradiction is now clear.
    print("\n--- Contradiction ---")
    print(f"Our initial choice (Assertion 1) was: f(x0) > {M_A}")
    print(f"Our derived conclusion is: f(x0) <= {M_A}")
    print("This is a logical contradiction.")
    print("\nTherefore, the initial assumption that the chain is positive recurrent must be false.")

prove_not_positive_recurrent()