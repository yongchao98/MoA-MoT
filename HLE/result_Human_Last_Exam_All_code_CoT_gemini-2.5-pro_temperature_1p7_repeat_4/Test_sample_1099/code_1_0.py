import math

def calculate_simulation_resources():
    """
    Calculates the minimal average resources needed to simulate the correlations
    of a singlet state in a CHSH-type scenario using a Local Hidden Variable model.
    """

    # --- 1. Communication Cost ---
    # According to Toner and Bacon (2003), to perfectly simulate the correlations
    # from projective measurements on a singlet state, 1 bit of communication is
    # necessary and sufficient on average.
    communication_cost_bits = 1

    # --- 2. Non-signaling PR-Box Cost ---
    # We can simulate the quantum correlations by mixing a maximally non-local
    # PR-Box with a classical (local) strategy.
    # Let p be the probability of using a PR-Box.
    # p * (PR-Box_score) + (1-p) * (Classical_limit) = Quantum_limit
    # p * 4 + (1-p) * 2 = 2 * sqrt(2)
    # Solving for p gives: p = sqrt(2) - 1

    tsirelson_bound = 2 * math.sqrt(2)
    classical_limit = 2.0
    pr_box_limit = 4.0

    # The equation to solve for p:
    # p = (tsirelson_bound - classical_limit) / (pr_box_limit - classical_limit)
    pr_box_cost_fraction = math.sqrt(2) - 1

    # --- 3. Shared Randomness Cost ---
    # The underlying local hidden variable in models that perfectly simulate the
    # singlet state (e.g., a shared random unit vector) requires infinite precision
    # to specify. Therefore, the amount of shared randomness is infinite.
    shared_randomness_cost = "infinite"

    print("Minimal Average Resources to Simulate Singlet State Correlations:")
    print("-" * 60)
    print(f"The equation for the PR-Box fraction (p) is derived from:")
    print(f"p * {pr_box_limit} + (1-p) * {classical_limit} = {tsirelson_bound:.4f}")
    print("\n--- Calculated Resource Costs ---")
    print(f"1. Average Communication Cost: {communication_cost_bits} bit")
    print(f"2. Average Non-signaling PR-Box Resource: {pr_box_cost_fraction:.4f}")
    print(f"3. Shared Randomness (LHV) Resource: {shared_randomness_cost} bits")
    print("-" * 60)

calculate_simulation_resources()

<<<1>>>