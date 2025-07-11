def calculate_agitation_medication():
    """
    Calculates and prints the recommended medication regimen for an agitated patient
    who has failed an initial dose of olanzapine.
    """
    # Initial medication that was already administered
    initial_olanzapine_mg = 5

    print(f"Clinical Situation: Patient remains agitated after an initial dose of {initial_olanzapine_mg}mg IM olanzapine.")
    print("The best next step is to use combination therapy for a synergistic effect.")
    print("-" * 30)

    # Based on the best choice (B), the following additional medications are given.
    additional_lorazepam_mg = 2
    additional_olanzapine_mg = 5
    
    # Calculate the new total dosages
    total_olanzapine_mg = initial_olanzapine_mg + additional_olanzapine_mg
    total_lorazepam_mg = 0 + additional_lorazepam_mg # Starting from 0 lorazepam

    print("Recommended Action (Choice B): Administer additional medication.")
    print(f"Dose to add: {additional_lorazepam_mg}mg IM lorazepam + {additional_olanzapine_mg}mg IM olanzapine")
    print("-" * 30)
    
    print("Total Dosage Calculation:")
    # Output each number in the final equation
    print(f"Total Olanzapine: {initial_olanzapine_mg}mg (initial) + {additional_olanzapine_mg}mg (additional) = {total_olanzapine_mg}mg")
    print(f"Total Lorazepam: 0mg (initial) + {additional_lorazepam_mg}mg (additional) = {total_lorazepam_mg}mg")
    print("-" * 30)
    
    print(f"The patient's final medication regimen will be: {total_olanzapine_mg}mg of olanzapine and {total_lorazepam_mg}mg of lorazepam.")

# Execute the function to display the answer
calculate_agitation_medication()
