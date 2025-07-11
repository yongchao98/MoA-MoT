def calculate_plant_fitness():
    """
    Models different insect behavioral patterns to determine which maximizes
    plant fitness, proxied by the number of flower interactions per hour.
    """
    total_time_seconds = 3600
    fitness_scores = {}
    
    print("Analyzing behavioral patterns based on their effect on plant fitness.")
    print("Plant fitness is estimated by the number of interactions (n(3)) per hour.\n")

    # --- Choice A: 4-3 >> 6-5 (Long interaction duration >> short feeding duration) ---
    # This implies long contact time per visit, but not necessarily many visits.
    t_investigate_A = 5  # seconds to investigate
    t_interact_A = 20    # seconds of a long interaction
    t_travel_A = 10      # seconds to travel to the next plant
    t_cycle_A = t_investigate_A + t_interact_A + t_travel_A
    # Each cycle contains one investigation and one interaction.
    n_cycles_A = total_time_seconds / t_cycle_A
    n_interactions_A = n_cycles_A * 1
    fitness_scores['A'] = n_interactions_A
    print(f"Pattern A (Long interaction duration):")
    print(f"  - Time per cycle (investigate + interact + travel) = {t_investigate_A} + {t_interact_A} + {t_travel_A} = {t_cycle_A}s")
    print(f"  - Number of interactions in 1hr = {total_time_seconds}s / {t_cycle_A}s/interaction = {n_interactions_A:.1f} interactions")
    print("-" * 30)

    # --- Choice B: 6-5 >> 4-3 ---
    # Logically impossible, as feeding (5->6) is a subset of interaction (3->4).
    fitness_scores['B'] = 0
    print("Pattern B (Feeding duration > Interaction duration):")
    print("  - Invalid. Feeding is a type of interaction; its duration cannot be greater.")
    print(f"  - Number of interactions in 1hr = {fitness_scores['B']} interactions")
    print("-" * 30)
    
    # --- Choice C: 4-3 >> 2-1 (Long interaction duration >> short investigation duration) ---
    # Similar to A, but with quicker decision-making.
    t_investigate_C = 2   # seconds for a short investigation
    t_interact_C = 20     # seconds for a long interaction
    t_travel_C = 10       # seconds to travel
    t_cycle_C = t_investigate_C + t_interact_C + t_travel_C
    n_cycles_C = total_time_seconds / t_cycle_C
    n_interactions_C = n_cycles_C * 1
    fitness_scores['C'] = n_interactions_C
    print(f"Pattern C (Long interaction, short investigation):")
    print(f"  - Time per cycle = {t_investigate_C} + {t_interact_C} + {t_travel_C} = {t_cycle_C}s")
    print(f"  - Number of interactions in 1hr = {total_time_seconds}s / {t_cycle_C}s/interaction = {n_interactions_C:.1f} interactions")
    print("-" * 30)
    
    # --- Choice D: n(5)/hour >> n(3)/hour ---
    # Logically impossible, you can't have more feeding starts than interaction starts.
    fitness_scores['D'] = 0
    print("Pattern D (Num feeding events > Num interaction events):")
    print("  - Invalid. A feeding event (5) must be part of an interaction event (3).")
    print(f"  - Number of interactions in 1hr = {fitness_scores['D']} interactions")
    print("-" * 30)

    # --- Choice E: n(1)/hour >> n(3)/hour (Many investigations, few interactions) ---
    # Insect investigates often but rarely makes contact. Bad for pollination.
    # Assume 1 interaction per 10 investigations.
    t_investigate_E = 5
    t_travel_E = 10
    # Time for 9 non-productive cycles (investigate + travel)
    time_non_productive = 9 * (t_investigate_E + t_travel_E)
    # Time for 1 productive cycle (investigate + interact + travel)
    time_productive = t_investigate_E + 10 + t_travel_E # Assume 10s interaction
    time_block_E = time_non_productive + time_productive
    interactions_per_block_E = 1
    n_blocks_E = total_time_seconds / time_block_E
    n_interactions_E = n_blocks_E * interactions_per_block_E
    fitness_scores['E'] = n_interactions_E
    print(f"Pattern E (Many investigations, few interactions):")
    print(f"  - Time per 10 investigations (with 1 interaction) = {time_block_E}s")
    print(f"  - Number of interactions in 1hr = ({total_time_seconds}s / {time_block_E}s) * {interactions_per_block_E} interaction = {n_interactions_E:.1f} interactions")
    print("-" * 30)
    
    # --- Choice F: n(3)/hour >> n(1)/hour (Many interactions per investigation) ---
    # This is the hallmark of an effective pollinator.
    # Once it finds a plant, it visits many flowers (many interactions).
    interactions_per_investigation_F = 10
    t_investigate_F = 5
    t_interact_single_F = 3  # Time for one short interaction (e.g., probing a flower)
    t_travel_F = 10
    # Time for one block: 1 investigation, 10 interactions, then travel
    t_block_F = t_investigate_F + (interactions_per_investigation_F * t_interact_single_F) + t_travel_F
    n_blocks_F = total_time_seconds / t_block_F
    n_interactions_F = n_blocks_F * interactions_per_investigation_F
    fitness_scores['F'] = n_interactions_F
    print(f"Pattern F (Many interactions per investigation):")
    print(f"  - Time for 1 investigation block with {interactions_per_investigation_F} interactions = {t_block_F}s")
    print(f"  - Number of interactions in 1hr = ({total_time_seconds}s / {t_block_F}s) * {interactions_per_investigation_F} interactions = {n_interactions_F:.1f} interactions")
    print("-" * 30)
    
    best_pattern = max(fitness_scores, key=fitness_scores.get)
    print(f"\nConclusion:")
    print(f"The behavior pattern '{best_pattern}' results in the highest number of interactions ({fitness_scores[best_pattern]:.1f}), and thus has the greatest positive effect on plant fitness.")

calculate_plant_fitness()