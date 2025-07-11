def demonstrate_violation():
    """
    This script demonstrates how critical-level views violate the
    Non-Elitism condition from Arrhenius's impossibility theorems.
    """

    # --- Step 1: Define the parameters for our scenario ---

    # A critical-level view states that the value of a population is the
    # sum of (utility - critical_level) for all individuals.
    # We choose a positive critical level.
    critical_level = 20

    # Non-Elitism considers an "elite" population with one person at a very high utility.
    elite_utility = 1000
    population_A = [elite_utility]

    # It also considers a population of many people at a low, but positive, utility.
    # CRUCIALLY, for this demonstration, we pick a low utility that is *below* the critical level.
    low_utility = 5  # Note: 5 > 0, but 5 < critical_level (20)

    # --- Step 2: Calculate the value of the 'elite' Population A ---

    # Value(A) = (utility_of_person_1 - critical_level)
    value_A = elite_utility - critical_level
    print(f"--- Scenario Setup ---")
    print(f"Critical Level (c): {critical_level}")
    print(f"Elite Population (A): One person with utility {elite_utility}")
    print(f"Low-Utility Population (B): 'N' people, each with utility {low_utility}\n")

    print(f"--- Value of Elite Population A ---")
    print(f"Value(A) = Elite Utility - Critical Level")
    print(f"Value(A) = {elite_utility} - {critical_level} = {value_A}\n")

    # --- Step 3: Show how the value of Population B changes with N ---

    print(f"--- Value of Low-Utility Population B ---")
    print("Non-Elitism requires that for some large number (N) of people in Population B,")
    print("Value(B) > Value(A). Let's test this.\n")

    # The value contribution of each person in Population B is negative
    value_per_person_B = low_utility - critical_level
    print(f"Value contributed by each person in B = {low_utility} - {critical_level} = {value_per_person_B}\n")

    # Let's check for increasing N
    for n in [1, 10, 100, 1000]:
        population_B = [low_utility] * n
        # Value(B) = N * (utility_of_each_person - critical_level)
        value_B = n * value_per_person_B
        is_better = value_B > value_A

        print(f"For N = {n}:")
        print(f"  Value(B) = {n} * ({low_utility} - {critical_level}) = {value_B}")
        print(f"  Is Value(B) > Value(A)? ({value_B} > {value_A}) -> {is_better}")

    print("\n--- Conclusion ---")
    print("As N increases, the value of Population B becomes more and more negative.")
    print("It can never become greater than the positive value of Population A.")
    print("Therefore, the Critical-Level View violates the Non-Elitism condition.")

# Run the demonstration
demonstrate_violation()
<<<C>>>