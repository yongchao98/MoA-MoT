import numpy as np

def simulate_mutation_rate_evolution():
    """
    This script simulates the evolution of the mutation rate in a population.
    It demonstrates that natural selection for an optimal fitness level
    maintains a relatively constant mutation rate.
    """
    # --- Simulation Parameters ---
    POPULATION_SIZE = 1000
    GENERATIONS = 50
    # This constant 'k' determines the strength of selection against high mutation rates.
    SELECTION_STRENGTH_K = 25.0
    # This constant determines how much the mutation rate itself can change per generation.
    MUTATION_RATE_EVOLVABILITY = 0.001

    # 1. Initialize a population with a wide, random distribution of mutation rates.
    # Each individual's primary trait for this simulation is its own mutation rate 'mu'.
    # We start with rates from 0.0 to 0.4.
    population_mu = np.random.rand(POPULATION_SIZE) * 0.4

    initial_avg_mu = np.mean(population_mu)
    print(f"Initial average mutation rate in the population: {initial_avg_mu:.5f}")

    # 2. Run the simulation over many generations.
    for gen in range(GENERATIONS):
        # 2a. Calculate fitness for each individual based on its mutation rate.
        # High mutation rates are penalized with lower fitness.
        fitness = np.exp(-SELECTION_STRENGTH_K * population_mu)

        # Normalize fitness values to get selection probabilities.
        total_fitness = np.sum(fitness)
        if total_fitness == 0:
            # If all fitnesses are zero (unlikely), all have an equal chance.
            probabilities = np.ones(POPULATION_SIZE) / POPULATION_SIZE
        else:
            probabilities = fitness / total_fitness

        # 2b. Select parents for the next generation based on fitness.
        # Individuals with higher fitness are more likely to be chosen.
        parent_indices = np.random.choice(
            a=np.arange(POPULATION_SIZE),
            size=POPULATION_SIZE,
            p=probabilities,
            replace=True
        )
        parent_mu = population_mu[parent_indices]

        # 2c. Create the next generation's mutation rates.
        # Offspring inherit the parent's mutation rate, plus a small random change.
        mutation_of_mu = np.random.normal(
            loc=0.0,
            scale=MUTATION_RATE_EVOLVABILITY,
            size=POPULATION_SIZE
        )
        population_mu = parent_mu + mutation_of_mu

        # Ensure mutation rates cannot become negative.
        population_mu[population_mu < 0] = 0

    final_avg_mu = np.mean(population_mu)
    print(f"Final average mutation rate after {GENERATIONS} generations: {final_avg_mu:.5f}")
    
    print("\nThe simulation shows selection favoring an optimal, low mutation rate.")
    print("This stabilizing pressure is driven by a fitness function.")

    # 3. Print the final equation that represents this selective pressure.
    print("\nFinal Equation representing the fitness-based selection:")
    print(f"Fitness = exp(-{SELECTION_STRENGTH_K} * Mutation_Rate)")

if __name__ == '__main__':
    simulate_mutation_rate_evolution()