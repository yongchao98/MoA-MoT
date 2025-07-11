def generate_counseling_recommendation():
    """
    Identifies a key drug interaction and generates a pharmacist's counseling recommendation.
    """
    # Patient's prescription medications and dosages
    drug_fluoxetine = "Fluoxetine"
    dose_fluoxetine = "20mg"
    drug_atorvastatin = "Atorvastatin"
    dose_atorvastatin = "20mg"
    drug_junel = "Junel Fe"
    dose_junel = "1.5/30mg"

    # Patient's OTC medication and its risky component
    otc_drug = "Excedrin"
    risky_component = "aspirin"

    # Safer alternative for pain relief
    safer_alternative = "acetaminophen (the main ingredient in Tylenol)"

    # Formulate the counseling message, including the specific drug dosages
    counseling_message = (
        f"Hi Allison, I have your prescriptions for {drug_atorvastatin} {dose_atorvastatin}, "
        f"{drug_junel} {dose_junel}, and {drug_fluoxetine} {dose_fluoxetine} ready.\n\n"
        f"I see from our chat that you took {otc_drug} for a headache. "
        f"I want to point out an important interaction with one of your new medicines.\n\n"
        f"Your {drug_fluoxetine} {dose_fluoxetine} can increase the risk of bleeding. "
        f"The {otc_drug} you took contains {risky_component}, which also increases this risk. "
        f"When taken together, they can significantly raise the chance of stomach bleeding.\n\n"
        f"For any future headaches, I recommend you choose a safer option like {safer_alternative}, "
        "as it does not have this interaction. Please feel free to ask any other questions."
    )

    print(counseling_message)

generate_counseling_recommendation()