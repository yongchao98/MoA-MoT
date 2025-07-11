def pharmacy_counseling():
    """
    Analyzes a patient's medications to identify interactions and generate a counseling recommendation.
    """
    # Patient's prescribed medications with their dosages
    prescriptions = {
        "Atorvastatin": 20,
        "Junel Fe": "1.5/30",
        "Fluoxetine": 20
    }
    # Patient's over-the-counter medication
    otc_medication = "Excedrin"

    # Known drug properties relevant to the interaction
    # Fluoxetine is an SSRI, which can increase bleeding risk.
    # Excedrin contains Aspirin, an NSAID, which also increases bleeding risk.
    ssri_drug = "Fluoxetine"
    nsaid_in_otc = "Aspirin"  # A key component of Excedrin

    recommendation = "There is a significant drug interaction risk. Please see the counseling recommendation below.\n"

    # Check if the patient is taking both an SSRI and an NSAID
    if ssri_drug in prescriptions:
        recommendation += f"\n"
        recommendation += "========================================\n"
        recommendation += "Pharmacist Counseling Recommendation\n"
        recommendation += "========================================\n\n"
        recommendation += f"Patient's medications considered:\n"
        recommendation += f"- Atorvastatin {prescriptions['Atorvastatin']}mg\n"
        recommendation += f"- Junel Fe {prescriptions['Junel Fe']}mg\n"
        recommendation += f"- Fluoxetine {prescriptions['Fluoxetine']}mg\n"
        recommendation += f"- Excedrin (contains {nsaid_in_otc})\n\n"

        recommendation += "--- IMPORTANT INTERACTION ---\n"
        recommendation += f"The {ssri_drug} you are taking should not be combined with Excedrin.\n\n"
        recommendation += f"1. The Problem: {ssri_drug} ({prescriptions[ssri_drug]}mg) can increase the risk of bleeding. Excedrin contains {nsaid_in_otc}, which also increases bleeding risk. Taking them together puts you at a much higher risk for stomach or intestinal bleeding.\n\n"
        recommendation += "2. The Recommendation: For headaches or other pain, please AVOID Excedrin, Aspirin, or Ibuprofen. A safer choice would be a medication containing only acetaminophen (such as Tylenol).\n\n"
        recommendation += "Please speak with the pharmacist if you have any questions."

    print(recommendation)

pharmacy_counseling()