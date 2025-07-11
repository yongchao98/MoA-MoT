def explain_violation():
    """
    This function explains and demonstrates which of Arrhenius's conditions
    is violated by critical-level views in population ethics.
    """
    print("--- Understanding the Concepts ---")
    print("1. Critical-Level View: A view in population ethics where adding a person to the world is only a net positive if their welfare 'w' is above a certain positive 'critical level' (c). The value contributed by that person is calculated as (w - c).")
    print("2. Non-Elitism Condition: This principle states that for any positive welfare level 'u', there must be some number of people 'n' such that a population of n people all at level u is better than a population of any number of people with extremely high welfare.")
    print("\n--- Demonstrating the Violation ---")

    # Define the parameters for our model
    critical_level_c = 20
    # We choose a welfare level 'u' that is positive, but below the critical level.
    welfare_level_u = 10
    # We create a comparison population 'B' with one person at a very high welfare level.
    population_b_welfare = 100

    print(f"Let's set a critical level c = {critical_level_c}.")
    print(f"Let's test the Non-Elitism condition with a welfare level u = {welfare_level_u}.")
    print(f"And let's create a comparison population B with one person at welfare = {population_b_welfare}.")

    # Calculate the value of population B
    value_of_b = population_b_welfare - critical_level_c
    print(f"\nThe value of population B is calculated as (welfare - c).")
    print(f"Value(B) = ({population_b_welfare} - {critical_level_c}) = {value_of_b}")

    print("\nNow, according to Non-Elitism, we should be able to find a number of people 'n' at level u=10 whose population is better than B.")

    # Calculate the value contributed by one person at level u
    value_per_person_at_u = welfare_level_u - critical_level_c
    print(f"\nLet's calculate the value of adding one person at welfare u = {welfare_level_u}.")
    print(f"Value contributed = ({welfare_level_u} - {critical_level_c}) = {value_per_person_at_u}")
    print("Since this value is negative, adding any number of people at this level will result in a negative total value.")

    # Let's try with a very large 'n'
    n = 1000
    total_value_of_n_population = n * value_per_person_at_u
    print(f"\nLet's test with a large n, for example n = {n}.")
    print(f"The total value of this population is n * (u - c).")
    print(f"Total Value = {n} * ({welfare_level_u} - {critical_level_c}) = {total_value_of_n_population}")

    print(f"\nIs the new population's value ({total_value_of_n_population}) greater than population B's value ({value_of_b})?")
    print(f"Result: {total_value_of_n_population > value_of_b}")

    print("\n--- Conclusion ---")
    print("No matter how large 'n' becomes, the total value will always be negative and thus can never be greater than the positive value of population B.")
    print("This demonstrates that critical-level views violate the 'Non-Elitism' condition.")

# Execute the demonstration
explain_violation()
<<<C>>>