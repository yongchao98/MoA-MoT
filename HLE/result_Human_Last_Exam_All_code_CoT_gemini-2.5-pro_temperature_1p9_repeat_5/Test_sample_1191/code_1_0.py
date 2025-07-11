import numpy as np

def simulate_fitness_tradeoff():
    """
    This simulation models the concept of an optimal mutation rate.
    Fitness is determined by a trade-off between two opposing costs:
    1.  The 'mutational load' (cost_of_mutations): The higher the mutation rate,
        the more deleterious mutations occur, reducing fitness.
    2.  The 'cost of fidelity': Achieving a very low mutation rate requires
        a significant metabolic investment in DNA repair machinery, which also
        reduces fitness.

    We find the mutation rate that maximizes overall fitness.
    """
    print("Simulating the effect of mutation rate on fitness...")

    # Constants for the model. These are for illustrative purposes.
    # c1 represents the severity of deleterious mutations.
    c1 = 50.0
    # c2 represents the metabolic cost of DNA fidelity/repair.
    c2 = 0.05

    # We will test mutation rates in this range
    possible_rates = np.linspace(0.01, 0.2, 200)

    # Variables to store the optimal point
    max_fitness = -float('inf')
    optimal_rate = 0
    cost_of_mutations_at_optimum = 0
    cost_of_fidelity_at_optimum = 0


    for rate in possible_rates:
        # Cost increases linearly with the mutation rate
        cost_of_mutations = c1 * rate
        # Cost increases as the rate gets closer to zero (hyperbolic)
        cost_of_fidelity = c2 / rate

        # Total fitness is a baseline value minus the two costs
        # A baseline of 10.0 is used to keep fitness positive in this model.
        current_fitness = 10.0 - cost_of_mutations - cost_of_fidelity

        if current_fitness > max_fitness:
            max_fitness = current_fitness
            optimal_rate = rate
            cost_of_mutations_at_optimum = cost_of_mutations
            cost_of_fidelity_at_optimum = cost_of_fidelity

    print("\nSimulation complete. The optimal mutation rate is found by balancing two costs.")
    print("This demonstrates 'Natural selection for fitness optimality'.")
    print("\n--- Optimal Point Found ---")
    print(f"Optimal Mutation Rate: {optimal_rate:.4f}")
    print("This rate maximizes fitness through a trade-off:")
    # Here we print the numbers that go into the final fitness calculation at the optimal point.
    print(f"Final Equation at Optimum: Fitness ({max_fitness:.4f}) = 10.0 - Deleterious_Cost ({cost_of_mutations_at_optimum:.4f}) - Fidelity_Cost ({cost_of_fidelity_at_optimum:.4f})")

simulate_fitness_tradeoff()