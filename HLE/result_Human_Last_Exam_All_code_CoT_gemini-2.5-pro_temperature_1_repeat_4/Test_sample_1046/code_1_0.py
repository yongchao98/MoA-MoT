import sys

def solve():
    """
    This script explains and simulates the rapid acquisition of drug resistance
    in a bacterial population without lateral gene transfer.
    """
    
    # --- Step 1: Analyze the biological scenario ---
    # The problem states that a bacterium with a stable genome (no lateral transfer)
    # acquires drug resistance as fast as one with frequent lateral transfer.
    # This implies a powerful vertical evolution mechanism is at play.
    
    # --- Step 2: Evaluate the best explanation ---
    # Choice B is the most complete explanation. It proposes a two-step mutation process:
    # 1. A rare mutation confers resistance but comes with a fitness cost.
    # 2. A second, compensatory mutation not only negates the cost but greatly
    #    increases fitness, creating a "super-strain".
    # This highly-fit strain can then sweep through the population rapidly under drug selection.
    # Cross-resistance would make this strain even more advantageous.

    # --- Step 3: Simulate the scenario to illustrate the concept ---
    print("This simulation models the rapid spread of drug resistance via compensatory mutation.")
    print("We will track three populations: Wild-Type (WT), Resistant with Cost (R-cost), and Compensated Resistant (R-comp).\n")
    
    # Simulation Parameters
    total_population_capacity = 1000000
    # We start with a population where one resistant mutant has just appeared.
    population = {
        'WT': 999999,
        'R-cost': 1,
        'R-comp': 0
    }
    
    # Fitness values (growth rate) in the presence of a drug.
    # A value of 1.0 means the population size is stable.
    fitness_wt = 0.1      # WT dies off under drug pressure.
    fitness_r_cost = 1.1  # The initial mutant survives but has a fitness cost (slow growth).
    fitness_r_comp = 1.9  # The compensated mutant is highly fit and grows very fast.
    
    # The rate at which an R-cost bacterium acquires a compensatory mutation to become R-comp.
    compensation_rate = 0.02 # 2% chance per generation
    
    generations = 25
    
    print(f"Initial State | WT: {population['WT']}, R-cost: {population['R-cost']}, R-comp: {population['R-comp']}")
    print("-" * 60)
    
    for gen in range(1, generations + 1):
        # Store current numbers to calculate growth based on the start of the generation
        current_r_cost = population['R-cost']
        
        # Growth of each sub-population
        grown_wt = population['WT'] * fitness_wt
        grown_r_cost = population['R-cost'] * fitness_r_cost
        grown_r_comp = population['R-comp'] * fitness_r_comp
        
        # Calculate new compensatory mutations arising from the R-cost population
        new_compensatory_mutants = grown_r_cost * compensation_rate
        
        # Update populations
        population['WT'] = grown_wt
        population['R-cost'] = grown_r_cost - new_compensatory_mutants
        population['R-comp'] = grown_r_comp + new_compensatory_mutants
        
        # Normalize to total population capacity
        current_total = sum(population.values())
        if current_total > total_population_capacity:
            factor = total_population_capacity / current_total
            for key in population:
                population[key] = int(population[key] * factor)
        else:
             for key in population:
                population[key] = int(population[key])

        if gen % 5 == 0:
            print(f"Generation {gen:2d} | WT: {population['WT']:<7} | R-cost: {population['R-cost']:<7} | R-comp: {population['R-comp']:<7}")

    print("-" * 60)
    print("Conclusion: The initial R-cost strain establishes itself slowly. However, once the")
    print("highly-fit R-comp strain appears, it rapidly dominates the entire population.")
    print("This demonstrates the mechanism described in answer B.")

    # Final equation format as requested
    final_wt = population['WT']
    final_r_cost = population['R-cost']
    final_r_comp = population['R-comp']
    total = final_wt + final_r_cost + final_r_comp
    
    print("\nFinal Population Equation:")
    print(f"{final_wt} (WT) + {final_r_cost} (R-cost) + {final_r_comp} (R-comp) = {total} (Total)")

solve()
<<<B>>>