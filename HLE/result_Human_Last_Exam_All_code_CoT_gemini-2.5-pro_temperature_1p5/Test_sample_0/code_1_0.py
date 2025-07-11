import math

def critical_level_value(n, w, c):
    """
    Calculates the value of a population using a critical-level view.
    Value = n * (w - c)
    n: number of people
    w: welfare level of each person
    c: the critical level
    """
    return n * (w - c)

def demonstrate_non_elitism_violation():
    """
    This function demonstrates how critical-level views violate Non-Elitism.
    Non-Elitism: For any w1 > w2 > 0, there exists a number n such that
                 a population of n people at w2 is at least as good as
                 one person at w1.
    Condition: value(n, w2, c) >= value(1, w1, c)
    """
    # --- Setup ---
    # Define a critical level, c.
    critical_level = 10.0

    # Define a high welfare level (w1) above the critical level.
    welfare_high = 20.0

    # Define a lower positive welfare level (w2) below the critical level.
    welfare_low = 5.0

    print("Demonstrating the violation of Non-Elitism by Critical-Level Views.")
    print(f"\nParameters Used:")
    print(f"  - Critical Level (c): {critical_level}")
    print(f"  - High Welfare Level (w1): {welfare_high}")
    print(f"  - Low Welfare Level (w2): {welfare_low}\n")

    # --- Value of the High-Welfare Individual ---
    # This value must be surpassed by a population at the low-welfare level.
    value_high_welfare_person = critical_level_value(1, welfare_high, critical_level)
    print("First, calculate the value of the one person at the high welfare level (w1).")
    print(f"Equation: 1 * (w1 - c)")
    print(f"Calculation: 1 * ({welfare_high} - {critical_level}) = {value_high_welfare_person}\n")

    # --- Test if the Condition Can Be Met ---
    print("Now, we check if any number of people 'n' at the low welfare level (w2)")
    print("can have a total value greater than or equal to this.")
    print(f"Condition to satisfy: n * (w2 - c) >= {value_high_welfare_person}")
    low_welfare_factor = welfare_low - critical_level
    print(f"Substituting values: n * ({welfare_low} - {critical_level}) >= {value_high_welfare_person}")
    print(f"Simplified: n * ({low_welfare_factor}) >= {value_high_welfare_person}\n")

    print("Since w2 is less than c, the term (w2 - c) is negative.")
    print("A positive number 'n' multiplied by a negative number will always be negative.")
    print("A negative result can never be greater than or equal to a positive one.")
    print("Let's test a few large values of 'n':\n")

    for n in [1, 100, 10000, 1000000]:
        value_low_welfare_pop = critical_level_value(n, welfare_low, critical_level)
        is_met = value_low_welfare_pop >= value_high_welfare_person
        print(f"--- For n = {n:,} ---")
        print(f"Equation: {n} * ({welfare_low} - {critical_level}) = {value_low_welfare_pop:,.1f}")
        print(f"Is {value_low_welfare_pop:,.1f} >= {value_high_welfare_person}? Result: {is_met}")

    print("\n--- CONCLUSION ---")
    print("No matter how large 'n' is, the total value for the low-welfare population remains")
    print("negative and thus can never meet or exceed the positive value of the single")
    print("high-welfare person. Therefore, critical-level views violate Non-Elitism.")

if __name__ == "__main__":
    demonstrate_non_elitism_violation()