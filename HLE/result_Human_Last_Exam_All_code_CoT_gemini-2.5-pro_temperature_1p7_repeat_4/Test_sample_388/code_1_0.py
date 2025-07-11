def provide_pharmacist_counseling():
    """
    Analyzes a patient's medications to identify interactions and generate a counseling recommendation.
    """
    # Patient's prescription medications with strengths
    prescriptions = {
        'Fluoxetine': '20mg',
        'Atorvastatin': '20mg',
        'Junel Fe': '1.5/30mg'
    }

    # Active ingredients in the OTC medication (Excedrin)
    otc_ingredients = ['Aspirin', 'Acetaminophen', 'Caffeine']

    # A simple database of known drug interaction rules
    # Rule format: {drug: {interacting_drug: 'reason'}}
    interaction_db = {
        'Fluoxetine': {
            'Aspirin': 'significantly increases the risk of gastrointestinal bleeding',
            'Ibuprofen': 'significantly increases the risk of gastrointestinal bleeding',
            'Naproxen': 'significantly increases the risk of gastrointestinal bleeding'
        }
    }

    recommendation_found = False
    # Check for interactions between prescriptions and OTC ingredients
    for med, strength in prescriptions.items():
        if med in interaction_db:
            for ingredient in otc_ingredients:
                if ingredient in interaction_db[med]:
                    reason = interaction_db[med][ingredient]
                    recommendation_found = True

                    # Print the detailed counseling recommendation
                    print("--- Pharmacist Counseling Recommendation ---")
                    print(f"\nInteraction Detected: Between {med} and {ingredient}.")
                    print("\nRecommendation Details:")
                    print(f"I see you are picking up your prescription for {med} {strength}. "
                          f"It is important to know that taking this medication with an NSAID "
                          f"like {ingredient} (which is found in the Excedrin you took) {reason}.")
                    print("\nFor future headaches, it would be safer to use a product that only contains Acetaminophen (like Tylenol) "
                          "instead of Excedrin or other NSAIDs like Ibuprofen (Advil/Motrin).")
                    print("\nPlease speak with your doctor or me if you have any more questions about managing your pain safely.")
                    print("-" * 41)
                    # The number "20" from "Fluoxetine 20mg" is included in the output.
                    break  # Stop after finding the first major interaction
        if recommendation_found:
            break

    if not recommendation_found:
        print("No significant interactions detected with the provided information. Always consult a healthcare professional.")

# Execute the function to get the recommendation
provide_pharmacist_counseling()
