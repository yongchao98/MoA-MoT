# Initial medication given
initial_olanzapine_mg = 5

# Proposed next step of medication (from the correct answer choice)
next_step_olanzapine_mg = 5
next_step_lorazepam_mg = 2

# Calculate the total dosage
total_olanzapine_mg = initial_olanzapine_mg + next_step_olanzapine_mg
total_lorazepam_mg = next_step_lorazepam_mg

# Print the breakdown of the final dosage
print("The patient initially received 5 mg of Olanzapine.")
print("The best next step is to administer an additional 5 mg of Olanzapine and 2 mg of Lorazepam.")
print("\nFinal Dosage Calculation:")
print(f"Total Olanzapine = {initial_olanzapine_mg} mg (initial) + {next_step_olanzapine_mg} mg (new) = {total_olanzapine_mg} mg")
print(f"Total Lorazepam = 0 mg (initial) + {next_step_lorazepam_mg} mg (new) = {total_lorazepam_mg} mg")