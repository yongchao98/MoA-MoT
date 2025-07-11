def demonstrate_violation():
    """
    This script demonstrates how critical-level views in population ethics
    violate the 'Non-Elitism' condition of Arrhenius's theorem.
    """

    # --- Setup ---
    # 1. Define a critical level greater than zero.
    critical_level = 10.0

    # 2. Define the 'Non-Elitism' condition to test.
    print("The 'Non-Elitism' condition states:")
    print("For any population of 'elites' with high welfare, it must be possible to imagine a different population with a vast number of people at a low (but positive) welfare level that is ethically better.")
    print("-" * 20)

    # --- Population A: The Elites ---
    # An elite population with one person enjoying very high welfare.
    elite_population_welfares = [100.0]
    elite_welfare = elite_population_welfares[0]

    # Calculate the value of the elite population.
    # Value = (welfare - critical_level)
    elite_value = (elite_welfare - critical_level)

    print("Step 1: Define an elite population A.")
    print(f"Population A has {len(elite_population_welfares)} person with welfare = {elite_welfare}")
    print(f"The critical level is {critical_level}.")
    print(f"Value of Population A = ({elite_welfare} - {critical_level}) = {elite_value}")
    print("-" * 20)

    # --- Population B: The Low-Welfare Population ---
    # An alternative population with many people at a low, but positive, welfare.
    # CRUCIAL: We pick a welfare level that is *below* the critical level.
    low_positive_welfare = 5.0
    num_people_in_b = 10000 # A very large number

    # Calculate the value contribution of a single person in this population.
    # Value per person = (welfare - critical_level)
    value_per_person_in_b = (low_positive_welfare - critical_level)

    # Calculate the total value of this large, low-welfare population.
    total_value_of_b = num_people_in_b * value_per_person_in_b

    print("Step 2: Try to create a better population B with low, positive welfare.")
    print(f"Let's choose a positive welfare level of {low_positive_welfare}, which is below the critical level of {critical_level}.")
    print(f"The value contribution of each person is ({low_positive_welfare} - {critical_level}) = {value_per_person_in_b}.")
    print(f"\nSince this value is negative, adding more people makes the population's total value worse.")
    print(f"Let's create Population B with {num_people_in_b} people.")
    print(f"Value of Population B = {num_people_in_b} * ({low_positive_welfare} - {critical_level}) = {total_value_of_b}")
    print("-" * 20)

    # --- Conclusion ---
    print("Step 3: Compare the populations.")
    print(f"Value of A (Elites): {elite_value}")
    print(f"Value of B (Many people): {total_value_of_b}")

    if total_value_of_b > elite_value:
        print("\nResult: Population B is better than Population A.")
    else:
        print("\nResult: Population B is NOT better than Population A.")
        print("\nConclusion: Because lives with positive welfare below the critical level contribute negatively, no number of such lives can outweigh the elite population. Therefore, critical-level views violate Non-Elitism.")

demonstrate_violation()