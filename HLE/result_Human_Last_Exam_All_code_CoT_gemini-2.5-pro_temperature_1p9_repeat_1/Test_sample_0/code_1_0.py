def demonstrate_violation():
    """
    This script demonstrates how critical-level views violate the Non-Elitism condition.
    """
    # 1. Define the critical level.
    # In critical-level views, this is a positive value.
    critical_level = 10

    # 2. Define Population A: a small population with very high welfare.
    pop_A_size = 1
    pop_A_welfare = 100 # Welfare is well above the critical level.

    # 3. Define Population B: a very large population with a welfare level
    # that is positive but below the critical level.
    pop_B_size = 1000000
    pop_B_welfare = 5 # Welfare is positive, but less than the critical level.

    # 4. Calculate the value of each population based on the critical-level formula.
    # The formula for total value is the sum of (welfare - critical_level) for all individuals.

    # Calculate value for Population A
    value_A = pop_A_size * (pop_A_welfare - critical_level)

    # Calculate value for Population B
    value_B = pop_B_size * (pop_B_welfare - critical_level)

    # 5. Print the results and explain the violation.
    print("--- Demonstrating a Violation in Population Ethics ---")
    print("The 'Non-Elitism' condition states that for any high-welfare population (A),")
    print("there should exist a larger, lower-positive-welfare population (B) that is better.\n")

    print(f"Let's set a positive critical level: {critical_level}\n")

    print("--- Population A (High Welfare) ---")
    print(f"Size: {pop_A_size}")
    print(f"Welfare per Person: {pop_A_welfare}")
    print("Value Calculation:")
    print(f"Equation: {pop_A_size} * ({pop_A_welfare} - {critical_level})")
    print(f"Total Value: {value_A}\n")

    print("--- Population B (Lower Positive Welfare) ---")
    print(f"Size: {pop_B_size}")
    print(f"Welfare per Person: {pop_B_welfare}")
    print("Value Calculation:")
    print(f"Equation: {pop_B_size} * ({pop_B_welfare} - {critical_level})")
    print(f"Total Value: {value_B}\n")

    print("--- Conclusion ---")
    print(f"Is Population B better than Population A? ({value_B} > {value_A}) -> {value_B > value_A}")
    print("Because Population B's welfare is below the critical level, each person adds negative value.")
    print("No matter how large Population B gets, its total value will remain negative")
    print("and can never surpass the positive value of Population A.")
    print("\nThis directly violates the condition of 'Non-Elitism'.")


if __name__ == '__main__':
    demonstrate_violation()
