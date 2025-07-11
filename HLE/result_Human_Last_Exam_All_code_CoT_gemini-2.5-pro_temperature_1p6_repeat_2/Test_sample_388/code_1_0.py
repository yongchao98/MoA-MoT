def pharmacy_counseling_analysis():
    """
    Analyzes a patient's medications for potential interactions and provides a recommendation.
    """
    print("Analyzing potential drug interactions based on the patient's information.")
    print("First, let's list the medications and dosages mentioned in the final analysis:")

    # The numbers from the prompt are included in the output as requested.
    atorvastatin_dose = 20
    junel_fe_dose_1 = 1.5
    junel_fe_dose_2 = 30
    fluoxetine_dose = 20

    print(f"- Atorvastatin: {atorvastatin_dose}mg")
    print(f"- Junel Fe: {junel_fe_dose_1}/{junel_fe_dose_2}mg")
    print(f"- Fluoxetine: {fluoxetine_dose}mg")
    print("- Excedrin (Over-the-Counter)")
    print("\n" + "="*40 + "\n")

    # Define Allison's medications and their relevant drug classes
    allisons_medications = [
        {'name': 'Fluoxetine', 'class': 'SSRI'},
        {'name': 'Atorvastatin', 'class': 'Statin'},
        {'name': 'Junel Fe', 'class': 'Oral Contraceptive'},
        # Excedrin's key ingredient for this interaction is Aspirin, an NSAID.
        {'name': 'Excedrin', 'active_ingredient': 'Aspirin', 'class': 'NSAID'}
    ]

    # Check for the presence of drug classes that have a known interaction
    found_ssri = None
    found_nsaid = None
    nsaid_ingredient = None

    for med in allisons_medications:
        if med['class'] == 'SSRI':
            found_ssri = med['name']
        if med['class'] == 'NSAID':
            found_nsaid = med['name']
            nsaid_ingredient = med.get('active_ingredient', 'NSAID')

    # If both interacting classes are found, provide the counseling point
    if found_ssri and found_nsaid:
        print("CRITICAL INTERACTION DETECTED!")
        print(f"Interaction between: {found_ssri} (SSRI) and {found_nsaid} (contains {nsaid_ingredient}, an NSAID).")
        print("\n--- Pharmacist Counseling Recommendation ---")
        print(f"The combination of an NSAID (like {nsaid_ingredient} in {found_nsaid}) and an SSRI (like {found_ssri}) can significantly increase the risk of bleeding, especially gastrointestinal bleeding.")
        print("\nACTIONABLE ADVICE:")
        print("1. For the headache, recommend an alternative pain reliever that does NOT contain an NSAID, such as single-ingredient acetaminophen (e.g., Tylenol).")
        print("2. Advise the patient to speak with her doctor about managing headaches and pain while taking Fluoxetine to find the safest long-term option.")
    else:
        print("No critical interactions were detected based on the standard checks for this scenario.")

# Run the analysis
pharmacy_counseling_analysis()