def provide_pharmacist_counseling():
    """
    Analyzes a patient's medications to identify key interactions and provide counseling.
    """
    # Patient's prescribed medications from the scenario
    prescriptions = {
        "Atorvastatin": "20mg",
        "Junel Fe": "1.5/30mg",
        "Fluoxetine": "20mg" # This is an SSRI
    }

    # Medications taken over-the-counter
    otc_medication = {
        "name": "Excedrin",
        "components": ["Acetaminophen", "Aspirin", "Caffeine"] # Aspirin is an NSAID
    }

    # Known drug interaction rule
    # SSRIs (like Fluoxetine) and NSAIDs (like Aspirin) can increase bleeding risk.
    has_ssri = "Fluoxetine" in prescriptions
    has_nsaid_in_otc = "Aspirin" in otc_medication["components"]

    print("Analyzing Allison's medications for potential interactions...")
    print("-" * 30)
    print(f"Prescriptions: {prescriptions}")
    print(f"OTC Medication: {otc_medication['name']} (contains {', '.join(otc_medication['components'])})")
    print("-" * 30)

    # If the specific interaction is found, provide the counseling recommendation.
    if has_ssri and has_nsaid_in_otc:
        fluoxetine_dose = prescriptions["Fluoxetine"]
        
        print("IMPORTANT COUNSELING RECOMMENDATION:")
        print(f"I see you are taking Fluoxetine {fluoxetine_dose}. It's important to know that taking Excedrin, which contains aspirin, at the same time can increase your risk of stomach or intestinal bleeding.")
        print("For future headaches, a safer option would be a product that only contains acetaminophen (like Tylenol) instead of a combination with aspirin.")
        print("Please be aware of symptoms of bleeding, such as unusual bruising, dark or tarry stools, or stomach pain. If you notice any of these, you should contact your doctor.")
    else:
        print("No critical interactions found based on the provided information.")

provide_pharmacist_counseling()