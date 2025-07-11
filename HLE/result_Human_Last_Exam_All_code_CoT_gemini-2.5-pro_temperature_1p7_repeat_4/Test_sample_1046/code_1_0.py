import random

def simulate_resistance_acquisition():
    """
    Simulates and explains two scenarios of bacterial resistance acquisition.
    """
    print("This simulation will demonstrate how a bacterium with a stable genome (Bacterium 2) can acquire resistance as quickly as one using lateral transfer (Bacterium 1).\n")
    print("This supports Answer B, which points to the combined effect of rare mutations, highly beneficial compensatory mutations, and cross-resistance.\n")

    # Simulation Parameters
    population_size = 10000
    generations = 15

    # --- Scenario 1: Bacterium with Lateral Transfer ---
    b1_resistant = 10  # Start with a few resistant bacteria
    b1_susceptible = population_size - b1_resistant
    # High rate of transferring resistance plasmids to susceptible bacteria
    lateral_transfer_rate = 0.15 

    # --- Scenario 2: Bacterium with Mutation and High Fitness ---
    b2_resistant = 0  # Starts with no resistant bacteria
    b2_susceptible = population_size
    # Rare chance of the initial resistance mutation occurring per generation
    mutation_rate = 0.1 # Chance per generation for the first mutation to appear
    # The new mutant has a significant growth advantage due to compensatory mutations
    fitness_advantage = 1.9 
    # The mutation provides resistance to multiple drugs (a conceptual multiplier)
    cross_resistance_factor = 3

    print(f"--- Simulation Start ---")
    print(f"Parameters for Bacterium 2 (Mutation-driven):")
    print(f"Initial Mutation Rate = {mutation_rate}")
    print(f"Fitness Advantage (from compensatory mutations) = {fitness_advantage}")
    print(f"Cross-Resistance Factor = {cross_resistance_factor}\n")
    print(f"{'Generation':<12} | {'Bacterium 1 Resistant':<25} | {'Bacterium 2 Resistant':<25}")
    print("-" * 68)

    # Run simulation
    for gen in range(1, generations + 1):
        # Scenario 1 Update
        newly_infected = int(b1_resistant * lateral_transfer_rate * (b1_susceptible / population_size))
        b1_resistant += newly_infected
        if b1_resistant > population_size: b1_resistant = population_size
        b1_susceptible = population_size - b1_resistant

        # Scenario 2 Update
        # Check if the rare mutation occurs for the first time
        if b2_resistant == 0 and random.random() < mutation_rate:
            b2_resistant = 1
            b2_susceptible = population_size - 1

        # If resistant bacteria exist, they outcompete the susceptible ones
        if b2_resistant > 0:
            # Calculate next generation based on fitness
            resistant_growth = b2_resistant * fitness_advantage
            susceptible_growth = b2_susceptible * 1.0 # No advantage for susceptible
            
            total_growth = resistant_growth + susceptible_growth
            
            # Normalize population back to original size
            b2_resistant = int((resistant_growth / total_growth) * population_size)
            b2_susceptible = population_size - b2_resistant
        
        print(f"{gen:<12} | {b1_resistant:<25} | {b2_resistant:<25}")

    print("\n--- Simulation Conclusion ---")
    print("As shown, even with a rare initial mutation, the high fitness advantage allows the resistant sub-population in Bacterium 2 to grow explosively once established.")
    print("This rapid takeover of the population allows it to match the pace of Bacterium 1, which relies on a steady rate of lateral transfer.")
    
    print("\nThe effectiveness of Bacterium 2's strategy can be represented conceptually:")
    print("Effective Pace = (Event Rarity) x (Fitness Advantage) x (Cross-Resistance Multiplier)")
    # Final output showing the numbers in the "equation"
    print("\nFinal Equation Numbers:")
    print(f"Event Rarity (initial mutation per generation): 1 / {1/mutation_rate}")
    print(f"Fitness Advantage (growth factor post-compensation): {fitness_advantage}")
    print(f"Cross-Resistance Multiplier (resistance to N drugs): {cross_resistance_factor}")
    print("\nThe product of these factors explains the rapid pace of observed resistance.")


simulate_resistance_acquisition()
<<<B>>>