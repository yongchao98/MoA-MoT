def calculate_agitation_medication():
    """
    Calculates and prints the total dosage for the recommended treatment plan.
    The patient has failed an initial dose of an antipsychotic. The best next step
    is to add a benzodiazepine and increase the antipsychotic to a full therapeutic dose.
    """
    initial_olanzapine_mg = 5
    additional_olanzapine_mg = 5
    lorazepam_mg = 2

    total_olanzapine_mg = initial_olanzapine_mg + additional_olanzapine_mg

    print("Clinical Recommendation: Add a benzodiazepine and increase the antipsychotic dose.")
    print(f"Initial Olanzapine Dose: {initial_olanzapine_mg} mg")
    print(f"Additional Olanzapine Dose: {additional_olanzapine_mg} mg")
    print(f"Lorazepam Dose to Add: {lorazepam_mg} mg")
    print("-" * 30)
    print("Final Plan:")
    print(f"Administer {lorazepam_mg} mg IM lorazepam + {additional_olanzapine_mg} mg olanzapine IM.")
    print(f"This results in a total olanzapine dose of {initial_olanzapine_mg} + {additional_olanzapine_mg} = {total_olanzapine_mg} mg.")

calculate_agitation_medication()