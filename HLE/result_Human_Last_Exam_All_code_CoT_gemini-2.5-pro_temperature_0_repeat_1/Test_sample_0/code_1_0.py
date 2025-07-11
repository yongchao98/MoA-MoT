def demonstrate_violation():
    """
    This script demonstrates which of Arrhenius's conditions is violated by
    critical-level views in population ethics.

    The condition we will test is Non-Elitism.
    - Non-Elitism states: There is some number 'n' of people with very low
      (but positive) welfare that is better than one person with very high welfare.
    - Critical-Level View states: The value of a population is the sum of
      (utility - critical_level) for each person.

    We will show that for a critical-level view, if the "low positive welfare"
    is below the critical level, Non-Elitism does not hold.
    """

    # 1. Define the parameters for our scenario
    critical_level = 10
    high_welfare = 100
    # Crucially, choose a low welfare level that is positive but BELOW the critical level
    low_welfare = 5
    # Let's test for a large number 'n' of low-welfare people
    number_of_low_welfare_people = 1000

    # 2. Calculate the value of each population according to the critical-level view
    # Value of one person with high welfare
    value_one_high_person = high_welfare - critical_level

    # Value of 'n' people with low welfare
    value_n_low_people = number_of_low_welfare_people * (low_welfare - critical_level)

    # 3. Print the explanation and results
    print("--- Testing Critical-Level Views against the Non-Elitism Condition ---")
    print(f"Parameters: Critical Level = {critical_level}, High Welfare = {high_welfare}, Low Welfare = {low_welfare}\n")

    print("Non-Elitism requires that for some number 'n', the value of 'n' people with low welfare")
    print("is greater than the value of one person with high welfare.\n")

    print("Calculating the value of one person at high welfare:")
    print(f"Value = (High Welfare - Critical Level)")
    print(f"Value = ({high_welfare} - {critical_level}) = {value_one_high_person}\n")

    print(f"Calculating the value of {number_of_low_welfare_people} people at low welfare:")
    print(f"Value = Number of People * (Low Welfare - Critical Level)")
    print(f"Value = {number_of_low_welfare_people} * ({low_welfare} - {critical_level}) = {value_n_low_people}\n")

    print("--- Conclusion ---")
    print(f"Is the value of the low-welfare population ({value_n_low_people}) greater than the value of the high-welfare person ({value_one_high_person})?")
    print(f"Result: {value_n_low_people > value_one_high_person}\n")

    print("Because the low welfare level is below the critical level, each person added contributes negative value.")
    print("Increasing the number of people only makes the total value more negative.")
    print("Therefore, a critical-level view violates the Non-Elitism condition.")

demonstrate_violation()
<<<C>>>