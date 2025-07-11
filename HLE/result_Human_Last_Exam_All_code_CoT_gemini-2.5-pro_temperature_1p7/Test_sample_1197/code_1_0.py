import random

def run_simulation(population_size, initial_frequency, selection_coefficient, generations):
    """
    Simulates the frequency of an allele over generations.
    A selection_coefficient of 0 represents pure genetic drift.
    A positive selection_coefficient represents adaptive evolution.
    """
    frequency = initial_frequency
    for gen in range(generations):
        # Stop if the allele is lost or fixed
        if frequency <= 0:
            return 0
        if frequency >= 1:
            return 1
            
        # 1. Selection step (if applicable)
        if selection_coefficient > 0:
            # Fitness of alleles: w_A (new allele) = 1 + s, w_a (old allele) = 1
            w_A = 1 + selection_coefficient
            w_a = 1
            
            # Calculate mean fitness of the population
            mean_fitness = frequency * w_A + (1 - frequency) * w_a
            
            # Update frequency based on selection
            frequency = (frequency * w_A) / mean_fitness
        
        # 2. Drift step (random sampling for the next generation)
        # Number of individuals is population_size. Diploid, so 2*N alleles.
        num_alleles = 2 * population_size
        
        # The number of new alleles in the next generation is a binomial draw
        # We can approximate this by drawing from a list
        alleles_in_next_gen = random.choices([1, 0], weights=[frequency, 1 - frequency], k=num_alleles)
        new_allele_count = sum(alleles_in_next_gen)
        
        # Update frequency based on the result of drift
        frequency = new_allele_count / num_alleles
        
    return frequency

def main():
    """Main function to run and print the simulation results."""
    # Simulation Parameters
    POP_SIZE = 500  # Size of the population
    INITIAL_FREQ = 0.01  # Starting frequency of the new allele (e.g., a single mutation)
    GENERATIONS = 300  # Number of generations to simulate
    NUM_RUNS = 10     # Number of independent simulations to run for each scenario

    print(f"--- Simulating Genetic Drift vs. Adaptive Evolution ---")
    print(f"Parameters: Population Size = {POP_SIZE}, Initial Frequency = {INITIAL_FREQ}, Generations = {GENERATIONS}")
    print("-" * 55)

    # --- Scenario 1: Pure Genetic Drift (Neutral Allele) ---
    print("\nScenario 1: Pure Genetic Drift (selection_coefficient = 0.0)")
    drift_fixations = 0
    for i in range(NUM_RUNS):
        final_freq = run_simulation(POP_SIZE, INITIAL_FREQ, 0.0, GENERATIONS)
        if final_freq == 1:
            drift_fixations += 1
        print(f"  Run {i+1}: Final frequency = {final_freq:.4f}")
    
    print(f"\nResult: Under pure drift, the allele reached fixation (frequency=1.0) in {drift_fixations} out of {NUM_RUNS} runs.")
    print("This outcome is random and fixation is rare for a new mutation.")
    print("-" * 55)

    # --- Scenario 2: Adaptive Evolution (Beneficial Allele) ---
    SELECTION_COEFF = 0.05  # 5% fitness advantage
    print(f"\nScenario 2: Adaptive Evolution (selection_coefficient = {SELECTION_COEFF})")
    selection_fixations = 0
    for i in range(NUM_RUNS):
        final_freq = run_simulation(POP_SIZE, INITIAL_FREQ, SELECTION_COEFF, GENERATIONS)
        if final_freq == 1:
            selection_fixations += 1
        print(f"  Run {i+1}: Final frequency = {final_freq:.4f}")

    print(f"\nResult: With a {SELECTION_COEFF*100}% fitness advantage, the allele reached fixation in {selection_fixations} out of {NUM_RUNS} runs.")
    print("This outcome shows how selection can predictably drive an allele's frequency, 'outweighing' the random effects of drift.")


if __name__ == "__main__":
    main()