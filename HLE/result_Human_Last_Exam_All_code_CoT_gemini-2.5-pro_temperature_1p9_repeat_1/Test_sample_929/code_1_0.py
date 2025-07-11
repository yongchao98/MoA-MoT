import sys

def solve():
    """
    This script models plant fitness based on insect behavior to determine
    which pattern is most beneficial for pollination.
    """

    # 1. Define fitness contribution rates (points per second)
    # Feeding is the most valuable action for pollination.
    R_feed = 10
    # Other forms of contact are less valuable but still useful.
    R_interact_non_feed = 1
    # Investigation (no contact) provides no fitness benefit.
    R_investigate = 0

    print("--- Plant Fitness Model ---")
    print(f"Fitness points per second for feeding (R_feed): {R_feed}")
    print(f"Fitness points per second for non-feeding interaction (R_interact_non_feed): {R_interact_non_feed}")
    print(f"Fitness points per second for investigation (R_investigate): {R_investigate}")
    print("-" * 40)
    print("Evaluating each behavioral pattern against a baseline...\n")

    results = {}

    def calculate_fitness(scenario_name, N_interact, D_interact, D_feed, D_invest, N_invest):
        """Calculates and prints the fitness score for a given scenario."""
        # Fitness from the feeding part of the interaction
        fitness_from_feed = N_interact * D_feed * R_feed
        # Fitness from the non-feeding part of the interaction
        D_interact_non_feed = D_interact - D_feed
        fitness_from_non_feed = N_interact * D_interact_non_feed * R_interact_non_feed
        # Total fitness from all interactions
        total_fitness = fitness_from_feed + fitness_from_non_feed
        
        # Total time is calculated for context but doesn't affect the fitness score itself
        total_time = (N_invest * D_invest) + (N_interact * D_interact)

        print(f"Scenario {scenario_name}:")
        print(f"  Parameters: N_interact={N_interact}, D_interact={D_interact}, D_feed={D_feed}, N_invest={N_invest}, D_invest={D_invest}")
        # The final equation as requested
        print(f"  Final Equation: (N_interact * D_feed * R_feed) + (N_interact * (D_interact - D_feed) * R_interact_non_feed)")
        print(f"  Calculation:   ({N_interact} * {D_feed} * {R_feed}) + ({N_interact} * ({D_interact} - {D_feed}) * {R_interact_non_feed}) = {total_fitness:.0f} points")
        print(f"  (Total behavior time: {total_time}s)")
        print("-" * 20)
        return total_fitness

    # --- Baseline Case (for comparison) ---
    base_N_interact = 50
    base_N_invest = 50
    base_D_interact = 20 # seconds per interaction
    base_D_feed = 10     # seconds of feeding per interaction
    base_D_invest = 10
    
    # --- A. 4-3 >> 6-5 (Interaction duration >> Feeding duration) ---
    # Interaction is long, but most of it is not feeding.
    A_D_interact = 100 
    A_D_feed = 10 
    results['A'] = calculate_fitness('A. 4-3 >> 6-5', base_N_interact, A_D_interact, A_D_feed, base_D_invest, base_N_invest)

    # --- B and D are logically impossible ---
    print("Scenarios B (6-5 >> 4-3) and D (n(5) >> n(3)) are logically impossible and are skipped.\n" + "-"*20)
    results['B'] = -1
    results['D'] = -1
    
    # --- C. 4-3 >> 2-1 (Interaction duration >> Investigation duration) ---
    # Long, committed interactions.
    C_D_interact = 100
    C_D_invest = 10
    C_D_feed = 50 # Assume feeding is proportional to interaction time
    results['C'] = calculate_fitness('C. 4-3 >> 2-1', base_N_interact, C_D_interact, C_D_feed, C_D_invest, base_N_invest)

    # --- E. n(1)/hr >> n(3)/hr (More investigations than interactions) ---
    # A "picky" insect that rarely makes contact.
    E_N_invest = 100
    E_N_interact = 10
    results['E'] = calculate_fitness('E. n(1) >> n(3)', E_N_interact, base_D_interact, base_D_feed, base_D_invest, E_N_invest)
    
    # --- F. n(3)/hr >> n(1)/hr (More interactions than investigations) ---
    # An "eager" insect that makes contact often.
    F_N_invest = 10
    F_N_interact = 100
    results['F'] = calculate_fitness('F. n(3) >> n(1)', F_N_interact, base_D_interact, base_D_feed, base_D_invest, F_N_invest)

    # --- Conclusion ---
    # Find the scenario with the maximum fitness score
    best_scenario = max(results, key=results.get)
    print("\n--- Conclusion ---")
    print(f"Scenario '{best_scenario}' results in the highest plant fitness score ({results[best_scenario]:.0f} points).")
    print("This pattern describes an efficient pollinator that frequently interacts with flowers, maximizing opportunities for pollination.")
    
solve()
<<<F>>>