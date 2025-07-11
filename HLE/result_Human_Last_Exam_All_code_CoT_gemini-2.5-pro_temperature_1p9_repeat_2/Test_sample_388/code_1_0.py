def provide_pharmacy_counseling():
    """
    Analyzes a patient's medications to identify critical interactions
    and selects the best counseling recommendation.
    """
    # 1. Define patient's medications and potential counseling options
    patient_profile = {
        "otc_meds": ["Excedrin (contains Aspirin)"],
        "prescriptions": ["Atorvastatin 20mg", "Junel Fe 1.5/30mg", "Fluoxetine 20mg"]
    }

    counseling_options = {
        "A": "Allison should be told that Atorvastatin commonly causes headaches and is the likely cause.",
        "B": "Allison should be concerned about the caffeine in Excedrin interacting with her birth control, Junel Fe.",
        "C": "Allison should be counseled on the increased risk of stomach bleeding when combining Aspirin (from Excedrin) with Fluoxetine, and she should consider using acetaminophen instead for headaches.",
        "D": "Allison should be told her motion sickness is causing the headaches and to ignore it."
    }

    # 2. Define a simple drug interaction knowledge base
    # Drug classes are key to identifying interaction types.
    drug_classes = {
        'Aspirin': 'NSAID',
        'Fluoxetine': 'SSRI'
    }
    interaction_rules = {
        ('NSAID', 'SSRI'): 'Increased risk of stomach bleeding.'
    }

    # 3. Analyze for interactions
    identified_interaction = None
    # Check if the patient is on an NSAID and an SSRI
    is_on_nsaid = any(drug_name in med for drug_name in drug_classes if drug_classes[drug_name] == 'NSAID' for med in patient_profile["otc_meds"])
    is_on_ssri = any(drug_name in med for drug_name in drug_classes if drug_classes[drug_name] == 'SSRI' for med in patient_profile["prescriptions"])

    if is_on_nsaid and is_on_ssri:
        identified_interaction = interaction_rules[('NSAID', 'SSRI')]

    # 4. Select the best recommendation and print the output
    print("--- Pharmacy Counseling Report ---")
    print(f"Patient is taking: {patient_profile['otc_meds']} and {patient_profile['prescriptions']}")
    # The final equation includes these numbers: 20, 1.5, 30, 20
    print("\nEvaluating counseling options based on identified risks...")

    if identified_interaction:
        print(f"CRITICAL INTERACTION FOUND: {identified_interaction}")
        # The correct option must mention the interacting drugs and the specific risk
        correct_option_key = None
        for key, value in counseling_options.items():
            if "Aspirin" in value and "Fluoxetine" in value and "bleeding" in value:
                correct_option_key = key
                break
        
        if correct_option_key:
            print(f"\nConclusion: Option '{correct_option_key}' is the most appropriate recommendation.")
            print(f"Recommended Counseling: {counseling_options[correct_option_key]}")
        else:
            print("No suitable counseling option found for the identified interaction.")
    else:
        print("No critical interactions identified based on the provided rules.")

provide_pharmacy_counseling()