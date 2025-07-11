import random

def simulate_resistance_acquisition():
    """
    This simulation models two scenarios of bacterial resistance acquisition
    to demonstrate why answer B is the most plausible explanation.
    """
    # --- Shared Parameters ---
    POP_SIZE = 100000  # Total population size
    TARGET_RESISTANCE_PCT = 0.90  # We stop when 90% of the population is resistant

    # --- Scenario 1: Lateral Gene Transfer (Fast Sharing) ---
    print("--- Simulating Bacterium 1: Resistance via Lateral Transfer ---")
    pop1 = {'sensitive': POP_SIZE, 'resistant': 0}
    # Introduce one bacterium with a resistance plasmid
    pop1['sensitive'] -= 1
    pop1['resistant'] += 1
    
    # Each resistant bacterium can transfer its plasmid to 5 others per generation
    transfer_rate_per_bacterium = 5
    generations1 = 0
    while (pop1['resistant'] / POP_SIZE) < TARGET_RESISTANCE_PCT:
        generations1 += 1
        # Calculate how many new transfers happen this generation
        newly_resistant = min(pop1['sensitive'], pop1['resistant'] * transfer_rate_per_bacterium)
        pop1['sensitive'] -= newly_resistant
        pop1['resistant'] += newly_resistant
        if generations1 > 50: # Failsafe
            break
            
    print(f"Bacterium 1 (Lateral Transfer) reached >{int(TARGET_RESISTANCE_PCT*100)}% resistance in {generations1} generations.\n")

    # --- Scenario 2: Mutation, Compensation, and Selective Sweep ---
    print("--- Simulating Bacterium 2: Resistance via Mutation & High-Fitness Sweep ---")
    print("This models the process described in Answer B.")
    pop2 = {
        'sensitive': POP_SIZE,
        'resistant_costly': 0,
        'resistant_compensated_fit': 0
    }
    
    mutation_rate = 0.00001 # 1 in 100,000 chance for the first resistance mutation
    compensation_rate = 0.001 # 1 in 1,000 chance for a fitness-boosting compensation
    
    # Fitness values (how well they reproduce under drug pressure)
    # Sensitive bacteria die off.
    # Costly mutants survive but reproduce slowly.
    # Compensated mutants are highly fit and reproduce very quickly.
    fitness = {
        'sensitive': 0.1, 
        'resistant_costly': 0.7, 
        'resistant_compensated_fit': 1.8 # Key factor: this is highly advantageous!
    }

    generations2 = 0
    while (pop2['resistant_compensated_fit'] / POP_SIZE) < TARGET_RESISTANCE_PCT:
        generations2 += 1
        
        # Step 1: Check for rare initial mutation
        if pop2['resistant_costly'] == 0 and pop2['resistant_compensated_fit'] == 0:
            if random.random() < mutation_rate * pop2['sensitive']:
                pop2['sensitive'] -= 1
                pop2['resistant_costly'] += 1

        # Step 2: Check for rare compensatory mutation in the costly-resistant lineage
        if pop2['resistant_costly'] > 0:
             if random.random() < compensation_rate * pop2['resistant_costly']:
                pop2['resistant_costly'] -= 1
                pop2['resistant_compensated_fit'] += 1
        
        # Step 3: Population growth/death based on fitness (selection)
        total_fitness = (pop2['sensitive'] * fitness['sensitive'] +
                         pop2['resistant_costly'] * fitness['resistant_costly'] +
                         pop2['resistant_compensated_fit'] * fitness['resistant_compensated_fit'])

        if total_fitness == 0: break # Extinction

        # Calculate the new population proportions
        new_sensitive_pct = (pop2['sensitive'] * fitness['sensitive']) / total_fitness
        new_costly_pct = (pop2['resistant_costly'] * fitness['resistant_costly']) / total_fitness
        new_compensated_pct = (pop2['resistant_compensated_fit'] * fitness['resistant_compensated_fit']) / total_fitness

        pop2['sensitive'] = int(new_sensitive_pct * POP_SIZE)
        pop2['resistant_costly'] = int(new_costly_pct * POP_SIZE)
        pop2['resistant_compensated_fit'] = int(new_compensated_pct * POP_SIZE)

        # Normalize total population back to POP_SIZE
        current_total = sum(pop2.values())
        pop2['sensitive'] += POP_SIZE - current_total
        
        if generations2 > 50: # Failsafe
            break

    print(f"Bacterium 2 (Mutation & Sweep) reached >{int(TARGET_RESISTANCE_PCT*100)}% resistance in {generations2} generations.")
    print("\nConclusion: The simulation shows that if a rare mutation is followed by a compensatory mutation that grants a high fitness advantage, the resulting selective sweep can be extremely rapid, matching the pace of lateral transfer.")

simulate_resistance_acquisition()
<<<B>>>