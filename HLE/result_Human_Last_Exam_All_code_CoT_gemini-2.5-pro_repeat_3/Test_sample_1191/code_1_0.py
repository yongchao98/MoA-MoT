import numpy as np

def illustrate_mutation_rate_optimality():
    """
    This simulation illustrates the principle that natural selection can
    drive the evolution of an optimal, constant genomic mutation rate.
    """
    print("Analyzing the effect of mutation rate on fitness...")
    print("The model assumes fitness is influenced by two opposing factors:")
    print("1. Genetic Load: High mutation rates increase harmful mutations, reducing fitness.")
    print("2. Adaptive Potential: A non-zero mutation rate is needed for adaptation, but its benefit saturates quickly.")
    print("-" * 50)

    # Define a range of possible genomic mutation rates (U)
    mutation_rates = np.linspace(0.01, 4, 400)

    # --- Fitness Calculation ---
    # 1. Fitness reduction from genetic load (e.g., fitness drops exponentially as U increases)
    fitness_from_load = np.exp(-mutation_rates)

    # 2. Small fitness benefit from adaptive potential (e.g., a small, saturating increase)
    fitness_from_adaptation = 0.05 * (1 - np.exp(-mutation_rates * 2))

    # Total fitness is the combination of these two effects
    total_fitness = fitness_from_load + fitness_from_adaptation

    # Find the optimal rate that maximizes fitness
    max_fitness_index = np.argmax(total_fitness)
    optimal_rate = mutation_rates[max_fitness_index]
    max_fitness = total_fitness[max_fitness_index]

    print(f"Simulation Result:")
    print(f"The mutation rate that maximizes fitness is approximately: {optimal_rate:.3f} mutations per genome per generation.")
    print(f"At this rate, the hypothetical fitness value is {max_fitness:.3f}.")
    print("\nThis demonstrates that natural selection would favor organisms with a mutation rate near this optimum.")
    print("Rates that are too high or too low result in lower fitness and would be selected against.")
    print("This selective pressure is the dominant factor that maintains an approximately constant genomic mutation rate across diverse species.")

# Run the illustrative simulation
if __name__ == "__main__":
    illustrate_mutation_rate_optimality()