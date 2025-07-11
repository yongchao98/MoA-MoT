def calculate_critical_level_value(population_welfares, critical_level):
    """
    Calculates the value of a population using a simple critical-level view.
    Value = Sum of (welfare - critical_level) for each person.
    """
    if not population_welfares:
        return 0
    return sum(w - critical_level for w in population_welfares)

# Step 1: Explain the principle being violated.
print("We will demonstrate how a critical-level view can violate the 'Weak Quality Addition' principle.")
print("Weak Quality Addition states: There is at least one scenario where adding a person with POSITIVE welfare makes a population better.")
print("-" * 70)

# Step 2: Set up a critical-level view that creates a violation.
# We do this by setting the critical level higher than the positive welfare level of the person we add.
critical_level_c = 20
welfare_of_new_person_u = 15

print(f"Let's define a critical-level view with a critical level c = {critical_level_c}.")
print(f"Now, let's consider adding one person with a POSITIVE welfare level u = {welfare_of_new_person_u}.")
print(f"Note that while the welfare u={welfare_of_new_person_u} is positive, it is below the critical level c={critical_level_c}.")
print("-" * 70)

# Step 3: Compare the population value before and after the addition.
# For simplicity, we start with an empty population (value = 0).
initial_population = []
initial_value = calculate_critical_level_value(initial_population, critical_level_c)

print(f"The value of the initial (empty) population is: {initial_value}")

# The new population has one person with welfare u.
new_population = [welfare_of_new_person_u]
new_value = calculate_critical_level_value(new_population, critical_level_c)

# Step 4: Show the calculation and the result.
print("\nThe value of the new population is calculated as (u - c).")
print(f"Final Value = ({welfare_of_new_person_u} - {critical_level_c}) = {new_value}")

print("\n--- Conclusion ---")
print(f"The value of the new population ({new_value}) is NOT greater than the initial population's value ({initial_value}).")
print("We have shown a case where adding a person with a positive welfare life did not make the population better.")
print("Because critical-level views permit this scenario, they violate the Weak Quality Addition principle.")

<<<E>>>