import random

def simulate_resistance_acquisition():
    """
    Simulates and compares two scenarios of resistance acquisition in bacteria.
    1. Population 1: Acquires resistance via Lateral Gene Transfer (LGT).
    2. Population 2: Acquires resistance via a rare mutation followed by a
       compensatory mutation that dramatically increases fitness and confers
       cross-resistance (as described in the correct answer choice).
    """

    # --- Simulation Parameters (The "Numbers in the Equation") ---
    print("--- Key Simulation Parameters ---")
    generations = 60
    print(f"Number of Generations = {generations}")
    
    # Population 1 (LGT) Parameters
    p_lgt_per_bacterium = 0.04 # 4% chance for a susceptible bacterium to receive a resistance gene each generation
    print(f"Probability of LGT per bacterium per generation = {p_lgt_per_bacterium}")

    # Population 2 (Mutation) Parameters
    p_initial_mutation = 0.0001  # Very rare event for resistance
    p_compensatory_mutation = 0.005 # A subsequent, less rare event
    print(f"Probability of initial resistance mutation = {p_initial_mutation}")
    print(f"Probability of compensatory mutation = {p_compensatory_mutation}")

    # Fitness values (relative growth rate)
    fitness_normal = 1.0
    fitness_costly_mutant = 0.8  # Initial resistance has a 20% fitness cost
    # The key: compensatory mutation leads to a massive fitness advantage
    fitness_super_mutant = 1.6
    print(f"Fitness of normal bacteria = {fitness_normal}")
    print(f"Fitness of initial resistant mutant = {fitness_costly_mutant}")
    print(f"Fitness of compensated 'super-mutant' = {fitness_super_mutant}\n")


    # --- Simulation for Population 1 (LGT) ---
    print("--- Population 1 (LGT) Results ---")
    pop1_proportions = {'susceptible': 1.0, 'resistant': 0.0}
    for gen in range(1, generations + 1):
        transfer = pop1_proportions['susceptible'] * p_lgt_per_bacterium
        pop1_proportions['susceptible'] -= transfer
        pop1_proportions['resistant'] += transfer
        if gen % 10 == 0 or gen == generations:
            print(f"Generation {gen}: {pop1_proportions['resistant']:.2%} of population is resistant.")


    # --- Simulation for Population 2 (Mutation) ---
    print("\n--- Population 2 (Mutation & Compensation) Results ---")
    pop2_proportions = {'susceptible': 1.0, 'costly_r': 0.0, 'super_r': 0.0}
    for gen in range(1, generations + 1):
        # Step 1: New mutations occur
        from_s_to_costly = pop2_proportions['susceptible'] * p_initial_mutation
        pop2_proportions['susceptible'] -= from_s_to_costly
        pop2_proportions['costly_r'] += from_s_to_costly

        from_costly_to_super = pop2_proportions['costly_r'] * p_compensatory_mutation
        pop2_proportions['costly_r'] -= from_costly_to_super
        pop2_proportions['super_r'] += from_costly_to_super

        # Step 2: Selection based on fitness (fitter types grow faster)
        total_fitness = (pop2_proportions['susceptible'] * fitness_normal +
                         pop2_proportions['costly_r'] * fitness_costly_mutant +
                         pop2_proportions['super_r'] * fitness_super_mutant)

        if total_fitness > 0:
            # Update proportions for the next generation based on who reproduces more
            pop2_proportions['susceptible'] = (pop2_proportions['susceptible'] * fitness_normal) / total_fitness
            pop2_proportions['costly_r'] = (pop2_proportions['costly_r'] * fitness_costly_mutant) / total_fitness
            pop2_proportions['super_r'] = (pop2_proportions['super_r'] * fitness_super_mutant) / total_fitness

        total_resistant = pop2_proportions['costly_r'] + pop2_proportions['super_r']
        if gen % 10 == 0 or gen == generations:
            print(f"Generation {gen}: {total_resistant:.2%} of population is resistant. (Dominated by super-mutant: {pop2_proportions['super_r']:.2%})")

if __name__ == '__main__':
    simulate_resistance_acquisition()