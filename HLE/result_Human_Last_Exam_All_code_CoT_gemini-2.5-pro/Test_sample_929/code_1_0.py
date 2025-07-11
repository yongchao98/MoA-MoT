def calculate_plant_fitness():
    """
    This function models the impact of different insect behavioral patterns on plant fitness.
    Fitness is primarily driven by pollination, which is most associated with feeding.

    Ethogram:
    1: investigation start, 2: investigation end
    3: interaction start,   4: interaction end
    5: feeding start,       6: feeding end

    We assign fitness weights to time spent in each activity:
    - Feeding: Highest value, as it directly leads to pollination.
    - Non-feeding Interaction: Low value, as pollination is possible but less likely.
    - Investigation: No value, as there is no contact.
    """
    # Fitness weights (points per unit of time)
    FITNESS_WEIGHT_FEEDING = 10.0
    FITNESS_WEIGHT_INTERACTION_NON_FEEDING = 1.0
    FITNESS_WEIGHT_INVESTIGATION = 0.0

    print("--- Evaluating Scenarios for Plant Fitness ---\n")

    # --- Scenario A: 4-3 >> 6-5 (Interaction duration >> Feeding duration) ---
    # Insect spends lots of time on the plant, but little time feeding.
    t_feeding_A = 10
    t_interaction_non_feeding_A = 90
    fitness_A = (t_feeding_A * FITNESS_WEIGHT_FEEDING) + (t_interaction_non_feeding_A * FITNESS_WEIGHT_INTERACTION_NON_FEEDING)
    print("Scenario A (4-3 >> 6-5):")
    print(f"  - Time spent feeding (6-5) = {t_feeding_A}")
    print(f"  - Time in non-feeding interaction = {t_interaction_non_feeding_A}")
    print(f"  - Equation: ({t_feeding_A} * {FITNESS_WEIGHT_FEEDING}) + ({t_interaction_non_feeding_A} * {FITNESS_WEIGHT_INTERACTION_NON_FEEDING}) = {fitness_A:.1f}")
    print(f"  - Resulting Fitness: {fitness_A:.1f}\n")

    # --- Scenario B: 6-5 >> 4-3 (Interpreted as Feeding duration is the majority of Interaction duration) ---
    # Insect is an efficient pollinator; most of its contact time is for feeding.
    t_feeding_B = 90
    t_interaction_non_feeding_B = 10
    fitness_B = (t_feeding_B * FITNESS_WEIGHT_FEEDING) + (t_interaction_non_feeding_B * FITNESS_WEIGHT_INTERACTION_NON_FEEDING)
    print("Scenario B (6-5 >> 4-3):")
    print(f"  - Time spent feeding (6-5) = {t_feeding_B}")
    print(f"  - Time in non-feeding interaction = {t_interaction_non_feeding_B}")
    print(f"  - Equation: ({t_feeding_B} * {FITNESS_WEIGHT_FEEDING}) + ({t_interaction_non_feeding_B} * {FITNESS_WEIGHT_INTERACTION_NON_FEEDING}) = {fitness_B:.1f}")
    print(f"  - Resulting Fitness: {fitness_B:.1f}\n")

    # --- Scenario C: 4-3 >> 2-1 (Interaction duration >> Investigation duration) ---
    # Insect quickly moves from investigating to interacting. We assume a balanced interaction.
    t_investigation_C = 10
    t_feeding_C = 50 # Assume interaction time is split 50/50 between feeding/non-feeding
    t_interaction_non_feeding_C = 50
    fitness_C = (t_feeding_C * FITNESS_WEIGHT_FEEDING) + (t_interaction_non_feeding_C * FITNESS_WEIGHT_INTERACTION_NON_FEEDING)
    print("Scenario C (4-3 >> 2-1):")
    print(f"  - Time spent investigating (2-1) = {t_investigation_C}")
    print(f"  - Time spent feeding (6-5) = {t_feeding_C}")
    print(f"  - Time in non-feeding interaction = {t_interaction_non_feeding_C}")
    print(f"  - Equation: ({t_feeding_C} * {FITNESS_WEIGHT_FEEDING}) + ({t_interaction_non_feeding_C} * {FITNESS_WEIGHT_INTERACTION_NON_FEEDING}) = {fitness_C:.1f}")
    print(f"  - Resulting Fitness: {fitness_C:.1f}\n")
    
    # --- Scenario E: n(1)/hour >> n(3)/hour (Investigation frequency >> Interaction frequency) ---
    # Many visitors, few make contact.
    # This is a frequency measure, but it implies low total interaction time.
    t_feeding_E = 5
    t_interaction_non_feeding_E = 5
    fitness_E = (t_feeding_E * FITNESS_WEIGHT_FEEDING) + (t_interaction_non_feeding_E * FITNESS_WEIGHT_INTERACTION_NON_FEEDING)
    print("Scenario E (n(1) >> n(3)):")
    print(f"  - Time spent feeding (6-5) = {t_feeding_E}")
    print(f"  - Time in non-feeding interaction = {t_interaction_non_feeding_E}")
    print(f"  - Equation: ({t_feeding_E} * {FITNESS_WEIGHT_FEEDING}) + ({t_interaction_non_feeding_E} * {FITNESS_WEIGHT_INTERACTION_NON_FEEDING}) = {fitness_E:.1f}")
    print(f"  - Resulting Fitness: {fitness_E:.1f}\n")


    print("--- Conclusion ---")
    print("Scenario B yields the highest fitness score because it maximizes the time spent in feeding, the behavior most crucial for pollination.")

calculate_plant_fitness()
<<<B>>>