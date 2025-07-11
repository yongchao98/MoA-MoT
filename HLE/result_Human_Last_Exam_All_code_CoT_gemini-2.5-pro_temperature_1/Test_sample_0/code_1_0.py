def check_non_elitism_violation():
    """
    This function demonstrates how Critical-Level Views violate
    the Non-Elitism condition from Arrhenius's impossibility theorem.
    """

    # --- Step 1: Define the parameters for our thought experiment ---
    # The critical level (c): a positive welfare level. Adding people below
    # this level is considered bad.
    critical_level_c = 10

    # The 'elite' welfare level (w): a high welfare level, above c.
    welfare_elite_w = 20

    # The 'lower' welfare level (x): a positive welfare level, but below c.
    welfare_lower_x = 5

    # The number of people in the 'elite' population.
    n_elite_population = 10

    print("--- Arrhenius Theorem Analysis: Critical-Level View vs. Non-Elitism ---")
    print(f"Setting up the scenario:")
    print(f"  - Critical Level (c) = {critical_level_c}")
    print(f"  - 'Elite' welfare level (w) = {welfare_elite_w}")
    print(f"  - 'Lower' welfare level (x) = {welfare_lower_x}")
    print(f"Note that we have chosen w > c > x, which is key to the violation.")
    print("-" * 60)

    # --- Step 2: Define the two populations as per the Non-Elitism condition ---
    print("Non-Elitism states that for a given Population A, there must exist a")
    print("sufficiently large Population B that is better.")
    print("\nLet's define our two populations:")
    print(f"  - Population A: {n_elite_population} people, each with welfare {welfare_elite_w}.")
    print(f"  - Population B: 'm' people, each with welfare {welfare_lower_x}.")
    print("-" * 60)

    # --- Step 3: Calculate the value of each population using the Critical-Level View ---
    print("According to the Critical-Level View, value is Sum(welfare - critical_level).")

    # Calculate the value of Population A
    value_A = n_elite_population * (welfare_elite_w - critical_level_c)
    print(f"\nCalculating value for Population A:")
    print(f"  Value(A) = {n_elite_population} * ({welfare_elite_w} - {critical_level_c})")
    print(f"  Value(A) = {n_elite_population} * ({welfare_elite_w - critical_level_c})")
    print(f"  Value(A) = {value_A}")

    # Calculate the per-person value contribution for Population B
    value_per_person_B = welfare_lower_x - critical_level_c
    print(f"\nCalculating value for Population B:")
    print(f"  The value added by each person in Population B is ({welfare_lower_x} - {critical_level_c}) = {value_per_person_B}.")
    print(f"  Therefore, Value(B) = m * {value_per_person_B}")
    print("-" * 60)

    # --- Step 4: Show the violation ---
    print("--- The Violation ---")
    print(f"Population A has a positive value ({value_A}).")
    print(f"Each person added to Population B contributes a negative value ({value_per_person_B}).")
    print("Therefore, no matter how large 'm' (the number of people in Population B) becomes,")
    print("the total value of Population B will always be negative.")

    print("\nA positive value can never be surpassed by a negative one.")
    print(f"So, Value(A) > Value(B) for any number of people m > 0.")

    print("\nConclusion: It is impossible to find a Population B that is better than Population A.")
    print("This directly violates the Non-Elitism condition.")


# Run the demonstration
check_non_elitism_violation()