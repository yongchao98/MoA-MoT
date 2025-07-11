import numpy as np

def model_fitness_vs_mutation_rate():
    """
    This function models the principle of an optimal mutation rate.
    It demonstrates how natural selection can maintain a stable rate by balancing
    the negative effects of a high rate (too many deleterious mutations) and
    a low rate (high metabolic cost of fidelity and lack of adaptability).
    """

    print("Modeling how natural selection finds an optimal mutation rate...\n")

    # --- Model Parameters ---
    # Represents the harm caused by an increasing number of deleterious mutations.
    # As mutation rate goes up, this cost increases linearly.
    cost_factor_deleterious = 50.0

    # Represents the metabolic cost of maintaining high-fidelity DNA replication/repair.
    # This cost is highest at very low mutation rates and decreases as the rate rises.
    cost_factor_fidelity = 0.01

    # We will test a range of possible mutation rates.
    mutation_rates = np.linspace(0.001, 0.05, 100)

    # Variables to store the best result found
    optimal_rate = 0
    max_fitness = -1

    print("Rate\t|\tFitness Score")
    print("----------------------------")

    # Iterate through each possible rate to find the one that maximizes fitness
    for rate in mutation_rates:
        # Cost from deleterious mutations (increases with rate)
        cost_deleterious = cost_factor_deleterious * rate

        # Cost from high fidelity (decreases with rate)
        cost_fidelity = cost_factor_fidelity / rate

        # Calculate a "fitness score". A higher score is better.
        # We start with a baseline and subtract the two costs.
        fitness_score = 1.0 - cost_deleterious - cost_fidelity

        # We only care about the peak; negative fitness isn't biologically real,
        # but the principle of finding the maximum score holds.
        if fitness_score > max_fitness:
            max_fitness = fitness_score
            optimal_rate = rate

    print(f"...\nSearch complete.\n")
    print("--- Conclusion ---")
    print(f"The mutation rate that results in the highest fitness is approximately: {optimal_rate:.4f}")
    print(f"This maximum fitness score is: {max_fitness:.4f}")
    print("\nThis model demonstrates that natural selection for fitness optimality (Choice A) can")
    print("act as a stabilizing force. It selects against rates that are too high (due to")
    print("deleterious mutations) and rates that are too low (due to metabolic costs and")
    print("lack of adaptation), thereby maintaining a relatively constant mutation rate.")

# Execute the model
model_fitness_vs_mutation_rate()