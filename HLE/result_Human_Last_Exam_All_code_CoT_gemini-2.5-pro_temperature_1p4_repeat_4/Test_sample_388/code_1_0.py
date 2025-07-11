def generate_counseling_recommendation():
    """
    Analyzes the patient's medications and generates a key counseling point.
    """
    
    # Medications involved in the primary interaction
    prescription_drug = "Fluoxetine"
    prescription_dose = 20
    otc_drug_component = "Aspirin (in Excedrin)"

    print("Pharmacist Counseling Recommendation:")
    print("-" * 35)

    print(f"The patient is taking {prescription_drug} {prescription_dose}mg and also taking {otc_drug_component} for a headache.")
    
    print("\nInteraction Alert:")
    print(f"Combining {prescription_drug} (an SSRI) with {otc_drug_component} (an NSAID) significantly increases the risk of stomach bleeding.")

    print("\nRecommendation:")
    print("1. It is recommended to avoid using Excedrin or other aspirin-containing products while taking Fluoxetine.")
    print("2. A safer alternative for the headache would be a medication containing only acetaminophen (e.g., Tylenol), as long as there are no contraindications.")
    print("3. The patient should be advised to watch for signs of gastrointestinal bleeding, such as dark or tarry stools, stomach pain, or vomiting blood, and to contact their doctor if these occur.")
    
    print("\nSecondary Point:")
    print("New or severe headaches can also be a side effect of hormonal contraceptives like Junel Fe 1.5/30mg. This should be discussed with her doctor as well.")


generate_counseling_recommendation()