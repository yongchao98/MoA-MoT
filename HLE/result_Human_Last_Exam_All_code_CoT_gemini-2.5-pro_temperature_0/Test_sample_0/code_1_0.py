def demonstrate_violation():
    """
    This function demonstrates how critical-level views violate the Non-Elitism condition
    of Arrhenius's impossibility theorem.
    """
    # 1. Define the parameters for our Critical-Level View example.
    # Let's set a positive critical level.
    critical_level = 10
    print(f"Let's assume a Critical-Level View where the critical level 'c' is {critical_level}.")
    print("According to this view, the value of a person with welfare 'w' is (w - c).\n")

    # 2. Define the parameters for the Non-Elitism test.
    # A single person with a very high quality of life.
    w_high = 1000
    # A group of people with a low, but positive, quality of life.
    # CRUCIALLY, we choose a welfare level that is below the critical level.
    w_low = 5
    # Let's test for a large number of people.
    n = 10000

    print("Now, let's test the 'Non-Elitism' condition.")
    print(f"Non-Elitism states that for a high welfare level (e.g., w_high = {w_high}),")
    print(f"there must be some number 'n' of people at a low positive welfare level (e.g., w_low = {w_low})")
    print("such that the group of 'n' people is considered better than the single high-welfare person.\n")

    # 3. Calculate the value of each population according to the critical-level view.
    value_high_welfare_person = w_high - critical_level
    value_low_welfare_group = n * (w_low - critical_level)

    print("Let's calculate the values:")
    print(f"Value of the single high-welfare person = w_high - c")
    print(f"Value = {w_high} - {critical_level} = {value_high_welfare_person}\n")

    print(f"Value of the group of n={n} low-welfare people = n * (w_low - c)")
    print(f"Value = {n} * ({w_low} - {critical_level}) = {n} * ({w_low - critical_level}) = {value_low_welfare_group}\n")

    # 4. Check if the Non-Elitism condition holds.
    # The condition is: value_low_welfare_group > value_high_welfare_person
    print("The Non-Elitism condition requires that:")
    print(f"Value of group > Value of single person")
    print(f"{value_low_welfare_group} > {value_high_welfare_person}\n")

    if value_low_welfare_group > value_high_welfare_person:
        print("Conclusion: The condition holds.")
    else:
        print("Conclusion: The condition is VIOLATED.")
        print(f"Because w_low ({w_low}) is less than the critical_level ({critical_level}), the term (w_low - c) is negative.")
        print("Therefore, as 'n' increases, the total value of the group becomes more negative and can never be greater than the positive value of the single high-welfare person.")

demonstrate_violation()
<<<C>>>