def calculate_agitation_protocol():
    """
    Calculates and displays the dosage for a severely agitated patient
    who has failed initial treatment.
    """
    initial_olanzapine_mg = 5
    
    # The chosen next step is option E
    additional_olanzapine_mg = 10
    additional_lorazepam_mg = 2
    
    # Calculate the total olanzapine administered
    total_olanzapine_mg = initial_olanzapine_mg + additional_olanzapine_mg
    
    print("Patient Scenario: A severely agitated and violent patient failed to respond to an initial dose.")
    print(f"Initial Dose Given: {initial_olanzapine_mg}mg IM olanzapine.")
    print("\nBest Next Step: Administer a combination of an antipsychotic and a benzodiazepine.")
    print("This provides a synergistic effect for rapid and safe sedation.")
    
    print("\nRecommended Additional Doses:")
    print(f"Olanzapine: {additional_olanzapine_mg}mg IM")
    print(f"Lorazepam: {additional_lorazepam_mg}mg IM")
    
    print("\nTotal Olanzapine Dosage Calculation:")
    print(f"{initial_olanzapine_mg}mg (initial) + {additional_olanzapine_mg}mg (additional) = {total_olanzapine_mg}mg total")

calculate_agitation_protocol()