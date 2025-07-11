import math

def calculate_mutant_fitness():
    """
    Models the overall fitness of a bacterial mutant that has acquired
    resistance, compensatory mutations, and cross-resistance.
    """
    # Start with the baseline fitness of a normal bacterium.
    base_fitness = 1.0

    # A resistance mutation often comes with a performance trade-off, reducing fitness.
    fitness_cost = 0.2  # Represents a 20% fitness cost.

    # A subsequent compensatory mutation can restore or even improve fitness.
    fitness_gain_from_compensation = 0.25 # Restores the 20% cost and adds a 5% bonus.

    # Cross-resistance means the mutation works against multiple drugs/stressors,
    # amplifying its selective advantage.
    cross_resistance_multiplier = 3.0 # The mutant is resistant to 3 different drugs.

    # Calculate the fitness of the new mutant strain
    # Net Fitness = (Base - Cost + Compensation) * Multiplier
    final_mutant_fitness_advantage = (base_fitness - fitness_cost + fitness_gain_from_compensation) * cross_resistance_multiplier

    print("This model demonstrates the powerful combination of factors leading to rapid resistance spread.")
    print("--------------------------------------------------------------------------------------")
    print(f"Base Fitness: {base_fitness}")
    print(f"Fitness Cost from initial resistance mutation: -{fitness_cost}")
    print(f"Fitness Gain from compensatory mutation: +{fitness_gain_from_compensation}")
    print(f"Advantage multiplier from cross-resistance: x{cross_resistance_multiplier}")
    print("--------------------------------------------------------------------------------------")
    print("The final equation for the mutant's selective advantage is:")
    print(f"({base_fitness} - {fitness_cost} + {fitness_gain_from_compensation}) * {cross_resistance_multiplier} = {final_mutant_fitness_advantage:.2f}")
    print("\nThis high fitness advantage allows the mutant strain to rapidly dominate the population.")

calculate_mutant_fitness()
<<<B>>>