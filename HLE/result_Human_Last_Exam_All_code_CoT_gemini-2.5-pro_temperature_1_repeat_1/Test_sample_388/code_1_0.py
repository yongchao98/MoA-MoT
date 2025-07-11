def provide_pharmacist_counseling():
    """
    Analyzes a patient's medications and symptoms to generate counseling recommendations.
    """

    # Patient's prescriptions
    prescriptions = {
        "Atorvastatin": "20mg",
        "Junel Fe": "1.5/30mg",
        "Fluoxetine": "20mg"
    }

    # Patient's OTC medication and its active ingredients
    otc_medication = "Excedrin"
    otc_components = ["Acetaminophen", "Aspirin", "Caffeine"]

    # Patient's reported symptom
    symptom = "headache"

    # --- Pharmacist's Analysis ---
    has_ssri = "Fluoxetine" in prescriptions
    has_nsaid = "Aspirin" in otc_components
    uses_oral_contraceptive = "Junel Fe" in prescriptions

    recommendations = []

    # Check for SSRI + NSAID interaction
    if has_ssri and has_nsaid:
        fluoxetine_dose = prescriptions["Fluoxetine"]
        recommendations.append(
            f"1. Be aware of an increased risk of stomach bleeding when taking {otc_medication} (which contains Aspirin) with your Fluoxetine {fluoxetine_dose}. "
            "Consider using an alternative like acetaminophen (Tylenol) for headaches in the future. Watch for signs of bleeding like unusual bruising or dark, tarry stools."
        )

    # Check for headache while on oral contraceptives
    if uses_oral_contraceptive and symptom == "headache":
        junel_fe_dose = prescriptions["Junel Fe"]
        recommendations.append(
            f"2. Since you are taking Junel Fe {junel_fe_dose}, it is important to discuss any new or worsening headaches with your doctor. "
            "They will want to ensure this is the best and safest contraceptive option for you."
        )
    
    # Print the full counseling recommendation
    print("Hello Allison, I have your prescriptions ready. I'd like to share a couple of important points about your medicines:")
    print("-" * 20)
    if recommendations:
        for rec in recommendations:
            print(rec)
    else:
        print("No specific interactions or immediate concerns noted with the information provided. Please take your medications as prescribed.")
    
    print("-" * 20)
    print("Let me know if you have any other questions.")

# Execute the function to print the counseling advice
provide_pharmacist_counseling()
