def demonstrate_violation():
    """
    This script demonstrates which condition of Arrhenius's theorem
    is violated by critical-level views in population ethics.
    """

    # 1. Define the Critical-Level View's parameters for our simulation.
    # The 'critical_level' is the welfare level a new life must exceed
    # to be considered a net positive addition to the world's value.
    critical_level = 10.0

    # 2. Define the 'Non-Elitism' principle's test case.
    # Non-Elitism states that for ANY positive welfare level, no matter how low,
    # there is SOME number of people whose existence is better than nothing.
    # We will test this with a low, positive welfare level that is BELOW the critical level.
    low_positive_welfare = 5.0

    print("--- Simulation: Testing Critical-Level Views against Non-Elitism ---")
    print(f"Setting a Critical Level (c) = {critical_level}")
    print(f"Using a low, positive welfare level (w) = {low_positive_welfare} for our test.")
    print("\nAccording to a critical-level view, the value each person adds is (w - c).")
    value_per_person = low_positive_welfare - critical_level
    print(f"Value added per person = {low_positive_welfare} - {critical_level} = {value_per_person}\n")

    print("Non-Elitism requires that we can find some number of people (n) to add")
    print("such that the total value is positive (better than an empty world with value 0).")
    print("Let's see if this is possible:\n")

    # 3. Demonstrate the conflict by testing with an increasing number of people (n).
    test_cases = [1, 10, 100, 1000000]
    for n in test_cases:
        total_value = n * value_per_person
        print(f"Adding {n} people:")
        # The final equation is n * (w - c)
        print(f"  Calculation: {n} * ({low_positive_welfare} - {critical_level}) = {total_value}")
        is_better_than_nothing = total_value > 0
        print(f"  Is the total value > 0? {is_better_than_nothing}\n")

    # 4. State the conclusion.
    print("--- Conclusion ---")
    print("As the simulation shows, because the welfare level (5.0) is below the critical level (10.0),")
    print("the value contributed by each person is negative.")
    print("Adding more people only makes the total value more negative; it never becomes positive.")
    print("\nThis directly violates the 'Non-Elitism' principle.")
    print("Therefore, critical-level views violate Non-Elitism.")

demonstrate_violation()
<<<C>>>