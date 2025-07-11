def generate_counseling_recommendation():
    """
    Analyzes the patient's medications and symptoms to generate a pharmacist's counseling recommendation.
    """
    patient_name = "Allison"
    prescriptions = {
        "Atorvastatin": "20mg",
        "Junel Fe": "1.5/30mg",
        "Fluoxetine": "20mg"
    }
    otc_med = "Excedrin (contains Aspirin)"
    symptom = "headache"

    # Main counseling points
    interaction = f"The Fluoxetine ({prescriptions['Fluoxetine']}) you are taking can increase the risk of bleeding when taken with the aspirin in Excedrin."
    recommendation_pain_relief = "For future headaches, a safer option would be a product containing only acetaminophen, like Tylenol."
    side_effect_warning = f"Also, new or worsening headaches can be a side effect of your birth control, Junel Fe ({prescriptions['Junel Fe']})."
    action_plan = "You should monitor this headache, and if it persists, worsens, or feels different from usual, it is important to contact your doctor."

    final_recommendation = (
        f"{interaction} {recommendation_pain_relief} {side_effect_warning} {action_plan}"
    )

    print("Generated Pharmacist Counseling Recommendation:")
    print(final_recommendation)


generate_counseling_recommendation()