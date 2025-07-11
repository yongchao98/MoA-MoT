def generate_counseling_recommendation():
    """
    Analyzes patient medication and generates a pharmacist's counseling recommendation.
    """
    # Patient's prescribed medications with dosages
    prescription_fluoxetine_name = "Fluoxetine"
    prescription_fluoxetine_dose = "20mg"
    prescription_fluoxetine_class = "SSRI"

    # Patient's over-the-counter medication
    otc_med_name = "Excedrin"
    otc_interacting_component = "aspirin"
    otc_component_class = "NSAID"

    # Identified risk and alternative
    primary_risk = "increased risk of stomach bleeding"
    safer_alternative = "acetaminophen (like Tylenol)"

    # Construct the recommendation message
    print("Pharmacist Counseling Recommendation:")
    print("=" * 35)
    
    # Using f-strings to embed variables and numbers into the output string
    recommendation_part_1 = f"Hello Allison. I see you are picking up your {prescription_fluoxetine_name} {prescription_fluoxetine_dose} today."
    recommendation_part_2 = f"It's important to know that taking {otc_med_name} with {prescription_fluoxetine_name} should be done with caution."
    recommendation_part_3 = f"The {otc_interacting_component} in {otc_med_name} is an {otc_component_class}, and when combined with an {prescription_fluoxetine_class} like {prescription_fluoxetine_name}, it can lead to an {primary_risk}."
    recommendation_part_4 = f"For future headaches, a safer choice would be a product containing only {safer_alternative}."
    recommendation_part_5 = f"Also, headaches can sometimes be a side effect of starting {prescription_fluoxetine_name}, so if they continue, please let your doctor know."

    print(recommendation_part_1)
    print(recommendation_part_2)
    print(recommendation_part_3)
    print(recommendation_part_4)
    print(recommendation_part_5)


# Execute the function to print the recommendation
generate_counseling_recommendation()