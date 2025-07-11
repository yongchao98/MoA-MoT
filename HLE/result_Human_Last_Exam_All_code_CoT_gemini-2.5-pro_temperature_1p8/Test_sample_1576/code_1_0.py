def calculate_agitation_medication():
    """
    Calculates and explains the appropriate next dose of medication for a
    severely agitated patient who failed initial treatment.
    """
    # Patient's initial treatment
    initial_medication = {'name': 'olanzapine', 'dose_mg': 5}

    # The chosen best next step involves a combination of two medications.
    # Choice B: 2mg IM lorazepam + 5mg olanzapine IM
    next_dose_olanzapine = 5
    next_dose_lorazepam = 2

    # Calculate the total dosage after the next step
    total_olanzapine = initial_medication['dose_mg'] + next_dose_olanzapine
    total_lorazepam = next_dose_lorazepam

    # Print the explanation and the final equation for total dosage
    print("Clinical Situation: A violent patient is inadequately sedated after an initial 5mg IM olanzapine.")
    print("The best next step is to administer a combination of an antipsychotic and a benzodiazepine.")
    print(f"\nRecommended next administration: {next_dose_olanzapine}mg olanzapine IM and {next_dose_lorazepam}mg lorazepam IM.")
    
    print("\n--- Final Dosage Calculation ---")
    print("This results in a total dosage of:")
    
    # Final equation for olanzapine
    print(f"Olanzapine: {initial_medication['dose_mg']}mg (initial) + {next_dose_olanzapine}mg (additional) = {total_olanzapine}mg total")
    
    # Final equation for lorazepam
    print(f"Lorazepam: 0mg (initial) + {next_dose_lorazepam}mg (additional) = {total_lorazepam}mg total")

# Execute the function to print the result
calculate_agitation_medication()