def generate_counseling_recommendation():
    """
    Analyzes Allison's medications and symptoms to provide a pharmacist's counseling recommendation.
    """

    # Prescription medications with dosages
    prescriptions = {
        "Atorvastatin": "20mg",
        "Junel Fe": "1.5/30mg",
        "Fluoxetine": "20mg"
    }

    # OTC medication used
    otc_medication = "Excedrin (contains Aspirin)"

    # Extracting the numbers for the final output string
    fluoxetine_dose = 20
    atorvastatin_dose = 20
    junel_fe_dose_1 = 1.5
    junel_fe_dose_2 = 30

    print("### Pharmacist Counseling Recommendation ###")
    print("\nBased on the medications being picked up and the information provided:")
    
    # Primary Counseling Point: Drug Interaction
    print("\n1. IMPORTANT - Medication Interaction Warning:")
    print(f"   Your prescription for Fluoxetine ({fluoxetine_dose}mg) can have a significant interaction with the aspirin in the Excedrin you took.")
    print("   - Both Fluoxetine and aspirin increase the risk of bleeding, especially in the stomach. Using them together raises this risk considerably.")
    print("   - RECOMMENDATION: For future headaches, please use a product that contains only acetaminophen (such as Tylenol) instead of Excedrin or other aspirin-containing products.")

    # Secondary Counseling Point: Headache Symptoms
    print("\n2. Physician Follow-up Recommended:")
    print(f"   You mentioned a new headache that started after a cruise and included an unusual symptom (smelling sea water).")
    print(f"   - New or worsening headaches can be a concern with certain medications, including your birth control (Junel Fe {junel_fe_dose_1}/{junel_fe_dose_2}mg).")
    print("   - RECOMMENDATION: It is important to mention these new symptoms to your doctor to ensure they are properly evaluated.")
    
    print("\nSummary of Medications:")
    print(f"- Atorvastatin {atorvastatin_dose}mg")
    print(f"- Junel Fe {junel_fe_dose_1}/{junel_fe_dose_2}mg")
    print(f"- Fluoxetine {fluoxetine_dose}mg")


# Run the function to print the recommendation
generate_counseling_recommendation()