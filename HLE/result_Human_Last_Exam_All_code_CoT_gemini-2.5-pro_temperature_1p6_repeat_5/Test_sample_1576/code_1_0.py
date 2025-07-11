# The patient has failed an initial intervention for severe agitation.
# This script calculates the total dosage of each medication after the recommended second intervention.

# Initial intervention
initial_olanzapine_mg = 5

# Second intervention based on the best choice (10mg IM olanzapine + 2mg IM lorazepam)
second_olanzapine_mg = 10
second_lorazepam_mg = 2

# Calculate the total amount of each drug administered
total_olanzapine_mg = initial_olanzapine_mg + second_olanzapine_mg
total_lorazepam_mg = second_lorazepam_mg # Lorazepam was not given initially

# Print the breakdown of the total dosages
print("Calculation for total Olanzapine (Zyprexa):")
print(f"{initial_olanzapine_mg} mg + {second_olanzapine_mg} mg = {total_olanzapine_mg} mg")

print("\nCalculation for total Lorazepam (Ativan):")
print(f"0 mg + {second_lorazepam_mg} mg = {total_lorazepam_mg} mg")
