def demonstrate_violation():
    """
    This script demonstrates how critical-level views violate the
    Non-Elitism principle from Arrhenius's impossibility theorems.
    """
    # 1. Define the parameters for our example
    critical_level = 10
    elite_welfare_u = 20
    non_elite_welfare_v = 5  # Note: u > critical_level > v
    elite_population_n = 1

    # 2. Define the valuation function for a critical-level view
    def get_population_value(size, welfare, crit_level):
        return size * (welfare - crit_level)

    # 3. Calculate the value of the 'elite' population A
    value_A = get_population_value(elite_population_n, elite_welfare_u, critical_level)

    print("--- Demonstrating the Violation of Non-Elitism ---")
    print(f"Critical Level (c): {critical_level}")
    print(f"Elite Welfare Level (u): {elite_welfare_u}")
    print(f"Non-Elite Welfare Level (v): {non_elite_welfare_v}")
    print("-" * 50)

    print("Step 1: Calculate the value of an 'elite' Population A.")
    print(f"Population A has {elite_population_n} person(s) at welfare level {elite_welfare_u}.")
    print("Value is calculated as: n * (u - c)")
    print(f"Value(A) = {elite_population_n} * ({elite_welfare_u} - {critical_level})")
    print(f"Value(A) = {value_A}\n")

    print("Step 2: Try to find a 'non-elite' Population B that is better than A.")
    print("Non-Elitism requires that for some size 'm', Population B is better than A.")
    print(f"Population B has 'm' people at welfare level {non_elite_welfare_v}.")
    print("Value is calculated as: m * (v - c)")

    # Let's test with a large 'm'
    non_elite_population_m = 1000
    value_B = get_population_value(non_elite_population_m, non_elite_welfare_v, critical_level)

    print(f"Let's test with m = {non_elite_population_m}:")
    print(f"Value(B) = {non_elite_population_m} * ({non_elite_welfare_v} - {critical_level})")
    print(f"Value(B) = {value_B}\n")

    print("Step 3: Check if Population B is better than Population A.")
    print("The condition to satisfy Non-Elitism is: Value(B) > Value(A)")
    print(f"The final equation we test is: {value_B} > {value_A}")

    if value_B > value_A:
        print("Result: This is TRUE. Non-Elitism is satisfied in this case.")
    else:
        print("Result: This is FALSE. Non-Elitism is violated.")

    print("\nConclusion: Because the welfare level 'v' is below the critical level 'c', its contribution is negative.")
    print("No matter how large 'm' is, the value of Population B will always be negative and can never be greater than the positive value of Population A.")
    print("This demonstrates that critical-level views violate the Non-Elitism principle.")

demonstrate_violation()