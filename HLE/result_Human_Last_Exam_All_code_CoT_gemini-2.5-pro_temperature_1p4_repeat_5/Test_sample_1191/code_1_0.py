import numpy as np

def analyze_mutation_rate_stability():
    """
    This model demonstrates how natural selection can maintain a stable mutation rate
    by balancing the fitness costs of too many mutations (mutational load)
    against the costs of preventing mutations (high fidelity).
    """

    print("--- Conceptual Model: Natural Selection on Mutation Rate ---")
    print("The hypothesis is that an optimal mutation rate emerges from a trade-off.")

    # 1. Cost of Mutational Load: Fitness decreases as the mutation rate (μ) increases,
    #    leading to more deleterious mutations. We model this as a direct cost.
    #    Equation: cost_load = C1 * μ
    C1 = 75.0  # A constant representing the negative impact of deleterious mutations.

    # 2. Cost of Fidelity: Maintaining a very low mutation rate is metabolically expensive
    #    (e.g., DNA repair machinery). We model this as an inverse cost.
    #    Equation: cost_fidelity = C2 / μ
    C2 = 0.003 # A constant representing the cost of replication/repair machinery.

    print("\nModel Parameters:")
    print(f"  - Load Cost Factor (C1): {C1}")
    print(f"  - Fidelity Cost Factor (C2): {C2}")

    # We search for the mutation rate (μ) that minimizes the total cost.
    # This represents the point of fitness optimality where natural selection would stabilize the rate.
    possible_mutation_rates = np.logspace(-4, -1, 500)
    min_total_cost = float('inf')
    optimal_rate = None
    optimal_load_cost = None
    optimal_fidelity_cost = None

    for rate in possible_mutation_rates:
        cost_load = C1 * rate
        cost_fidelity = C2 / rate
        total_cost = cost_load + cost_fidelity

        if total_cost < min_total_cost:
            min_total_cost = total_cost
            optimal_rate = rate
            optimal_load_cost = cost_load
            optimal_fidelity_cost = cost_fidelity

    print(f"\nSearching for the optimal mutation rate (μ) that minimizes total cost...")
    print(f"The predicted optimal mutation rate is approximately: {optimal_rate:.6f}")

    print("\nThis optimal rate is found at the point where:")
    print(f"  Total Fitness Cost = Cost of Load + Cost of Fidelity")
    # Here we output the numbers in the final calculation for the optimal point.
    print(f"  {min_total_cost:.4f} = {optimal_load_cost:.4f} + {optimal_fidelity_cost:.4f}")

    print("\nThis model illustrates 'Natural selection for fitness optimality' (Choice A),")
    print("which is the dominant factor in maintaining a relatively constant genomic mutation rate.")

analyze_mutation_rate_stability()