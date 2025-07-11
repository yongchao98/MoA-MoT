def calculate_critical_level_value(welfare_levels, critical_level):
    """Calculates the value of a population based on a critical-level view."""
    total_value = 0
    for welfare in welfare_levels:
        total_value += (welfare - critical_level)
    return total_value

# --- Step 1: Define the parameters for our example ---
critical_level = 10
high_welfare = 100
low_welfare = 5  # Note: This is > 0, but < critical_level

# --- Step 2: Define two populations to compare ---
# Population A: A small "elite" population
pop_A_welfare = [high_welfare]
pop_A_size = len(pop_A_welfare)

# Population B: A very large population with low welfare
# According to Non-Elitism, if we make this population big enough,
# it should be better than Population A.
pop_B_size = 1000000
pop_B_welfare = [low_welfare] * pop_B_size

# --- Step 3: Calculate the value of each population using the Critical-Level View ---
value_A = calculate_critical_level_value(pop_A_welfare, critical_level)
value_B = calculate_critical_level_value(pop_B_welfare, critical_level)

# --- Step 4: Print the results and the explanation ---
print("Demonstrating how Critical-Level Views violate Non-Elitism:\n")
print(f"Let's assume a Critical Level of: {critical_level}\n")

print("--- Population A (The 'Elite' Population) ---")
print(f"It has {pop_A_size} person with a welfare level of {high_welfare}.")
print(f"The value equation is: {pop_A_size} * ({high_welfare} - {critical_level})")
print(f"Calculated Value of Population A: {value_A}\n")

print("--- Population B (The Large, Lower-Welfare Population) ---")
print(f"It has {pop_B_size} people, each with a welfare level of {low_welfare}.")
print("Non-Elitism states that for a large enough population, this should be better than Population A.")
print(f"The value equation is: {pop_B_size} * ({low_welfare} - {critical_level})")
print(f"Calculated Value of Population B: {value_B}\n")

print("--- Conclusion ---")
print(f"The value of Population A is {value_A}, while the value of Population B is {value_B}.")
print("Because Population B's value is negative, it can never be better than Population A's, no matter how large it gets.")
print("This is a direct violation of the Non-Elitism condition.")
print("\nTherefore, the correct answer is C.")
