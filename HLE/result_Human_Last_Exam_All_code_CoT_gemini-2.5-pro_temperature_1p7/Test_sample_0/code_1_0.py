def demonstrate_violation():
    """
    This script demonstrates how Critical-Level Views violate the
    Non-Elitism condition of Arrhenius's sixth impossibility theorem.
    """

    # --- Step 1: Define parameters for our model ---
    # The critical level: any life below this has negative value.
    critical_level = 20
    # A low but positive welfare level, chosen to be below the critical level.
    low_positive_welfare = 5
    # A very high welfare level for one person.
    high_welfare = 1000
    # A large number of people for the Non-Elitism test.
    large_number_of_people = 100000

    # --- Step 2: Define the utility calculation for critical-level views ---
    def critical_level_value(welfare):
        """Calculates the value of a single person's life."""
        return welfare - critical_level

    # --- Step 3: Analyze the two populations ---

    print("--- The Setup ---")
    print(f"Critical Level (c): {critical_level}")
    print(f"High Welfare level (h): {high_welfare}")
    print(f"Low Positive Welfare level (w): {low_positive_welfare} (note: w > 0, but w < c)")
    print("-" * 20)

    # Population A: One person with very high welfare
    print("Population A: One person with high welfare.")
    value_A = critical_level_value(high_welfare)
    print(f"The value equation is: h - c")
    print(f"The final equation is: {high_welfare} - {critical_level} = {value_A}")
    print("-" * 20)

    # Population B: Many people with low, positive welfare
    print(f"Population B: {large_number_of_people:,} people with low positive welfare.")
    # First, calculate the value contributed by a single person in this population.
    value_per_person_B = critical_level_value(low_positive_welfare)
    print(f"The value equation for each person is: w - c")
    print(f"The final equation for each person is: {low_positive_welfare} - {critical_level} = {value_per_person_B}")

    # Now, calculate the total value of Population B.
    total_value_B = large_number_of_people * value_per_person_B
    print(f"The total value equation for Population B is: number_of_people * (w - c)")
    print(f"The final equation for total value is: {large_number_of_people} * ({low_positive_welfare} - {critical_level}) = {total_value_B}")
    print("-" * 20)

    # --- Step 4: Conclusion ---
    print("--- Conclusion ---")
    print("Non-Elitism requires that for some large number of people, Population B should be better than Population A.")
    print(f"Is Population B's value ({total_value_B}) greater than Population A's value ({value_A})?")
    print(f"Result: {total_value_B > value_A}")
    print("\nBecause the welfare of people in Population B is below the critical level, each person adds negative value.")
    print("Therefore, adding more people only makes the total value worse. Population A will always be considered better.")
    print("This directly violates the Non-Elitism condition.")


demonstrate_violation()
<<<C>>>