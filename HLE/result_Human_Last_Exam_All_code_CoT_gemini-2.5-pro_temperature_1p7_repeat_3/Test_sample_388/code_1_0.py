def pharmacist_counseling():
    """
    This function generates a counseling recommendation based on a patient's medications.
    """
    patient_name = "Allison"
    otc_medication = "Excedrin"
    prescription_med_1 = "Atorvastatin 20mg"
    prescription_med_2 = "Junel Fe 1.5/30mg"
    prescription_med_3 = "Fluoxetine 20mg"
    
    # Key components for interaction check
    fluoxetine = "Fluoxetine (an SSRI)"
    excedrin_component = "aspirin (an NSAID found in Excedrin)"
    interaction_risk = "increased risk of stomach bleeding"
    alternative_pain_reliever = "a product with just acetaminophen, like Tylenol"

    recommendation = (
        f"Hello {patient_name}, I am counseling you on your prescriptions today. "
        f"I see you are picking up your {prescription_med_3}. You mentioned taking {otc_medication} for your headache.\n\n"
        f"It is important for you to know that taking {fluoxetine} together with {excedrin_component} "
        f"can significantly increase the {interaction_risk}. \n\n"
        f"For any future headaches, a safer choice for you would be {alternative_pain_reliever}. "
        "Please let me know if you have any questions about this."
    )
    
    print(recommendation)

pharmacist_counseling()