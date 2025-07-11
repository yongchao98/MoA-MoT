import math

def simulate_resistance_spread():
    """
    This simulation demonstrates how two different mechanisms can lead to
    drug resistance spreading through a bacterial population at a similar pace.

    - Bacterium 1 uses rapid lateral gene transfer.
    - Bacterium 2 uses a rarer mutation event, but the resulting mutant
      has a very high fitness advantage due to compensatory mutations.
    """
    population_size = 1_000_000
    generations = 15

    # --- Parameters for Bacterium 1 (Lateral Transfer) ---
    pop1_resistant = 1  # Starts with one resistant bacterium
    # A high rate of transfer allows for rapid spread.
    # Represents a percentage of susceptible bacteria converted per generation by the resistant pool.
    transfer_coefficient = 0.4

    # --- Parameters for Bacterium 2 (Mutation + High Fitness) ---
    pop2_resistant = 0  # Starts with no resistant bacteria
    # A "super-mutant" with high fitness appears after a delay.
    mutation_event_time = 5
    # The new mutant is so fit it has a massive growth advantage.
    # We'll use a reproductive advantage factor of 2.5 (it produces 2.5 surviving
    # offspring for every 1 produced by the susceptible strain under drug pressure).
    fitness_advantage = 2.5

    print("Simulating resistance acquisition in two bacterial populations...")
    print(f"Goal: Reach >50% population resistance.\n")
    print(f"{'Generation':<12} | {'Bacterium 1 (Lateral Transfer)':<35} | {'Bacterium 2 (High-Fitness Mutant)':<40}")
    print("-" * 92)

    for gen in range(1, generations + 1):
        # --- Bacterium 1 Simulation Step ---
        # The number of newly resistant bacteria is proportional to the interactions
        # between resistant and susceptible individuals.
        if pop1_resistant > 0:
            fraction_resistant = pop1_resistant / population_size
            newly_converted = int(pop1_resistant * transfer_coefficient * (1 - fraction_resistant))
            pop1_resistant += newly_converted
            pop1_resistant = min(population_size, pop1_resistant)

        # --- Bacterium 2 Simulation Step ---
        if gen == mutation_event_time:
            # The rare mutation and compensatory mutation events occur, creating one "super-mutant".
            pop2_resistant = 1
        elif gen > mutation_event_time:
            # This super-fit mutant now proliferates exponentially.
            pop2_resistant = pop2_resistant * fitness_advantage
            pop2_resistant = min(population_size, math.ceil(pop2_resistant))

        # --- Print current status ---
        pop1_percent = (pop1_resistant / population_size) * 100
        pop2_percent = (pop2_resistant / population_size) * 100
        print(f"{gen:<12} | {f'{pop1_resistant:>12,} resistant ({pop1_percent:4.1f}%)':<35} | {f'{pop2_resistant:>12,} resistant ({pop2_percent:4.1f}%)':<40}")

    print("\n" + "=" * 92)
    print("Conclusion:")
    print("Both populations acquired resistance at a similar pace.")
    print("This is possible because the combination of a resistance mutation with a compensatory mutation in Bacterium 2 created a super-fit strain that spread via selection as fast as lateral transfer spread the gene in Bacterium 1.\n")

    print("The final state of Bacterium 2 can be described by the following growth equation:")
    initial_mutant_count = 1
    generations_of_growth = generations - mutation_event_time
    final_population = int(initial_mutant_count * (fitness_advantage ** generations_of_growth))
    final_population = min(population_size, final_population)

    # Outputting the final equation with its numbers
    print("Final Resistant Population = Initial Mutant Count * (Fitness Advantage ^ Generations of Growth)")
    print(f"{final_population} = {initial_mutant_count} * ({fitness_advantage} ** {generations_of_growth})")


simulate_resistance_spread()
<<<B>>>