def provide_counseling():
    """
    Analyzes a patient's medications to identify potential interactions and
    generates a pharmacist's counseling recommendation.
    """
    # Define patient's medications and their relevant classes
    prescriptions = {
        'Fluoxetine': {'class': 'SSRI'},
        'Atorvastatin': {'class': 'Statin'},
        'Junel Fe': {'class': 'Oral Contraceptive'}
    }

    # Define the OTC medication and its active components
    otc_drug = 'Excedrin'
    otc_components = {
        'Aspirin': {'class': 'NSAID'},
        'Acetaminophen': {'class': 'Analgesic'},
        'Caffeine': {'class': 'Stimulant'}
    }

    # Define a simple database of known clinically significant interactions
    # A frozenset is used for the key so the order of drug classes does not matter.
    interaction_db = {
        frozenset(['SSRI', 'NSAID']): (
            "Combining an SSRI with an NSAID can significantly increase the "
            "risk of stomach and intestinal bleeding."
        )
    }

    counseling_recommendation = ""

    # Check for interactions between prescribed and OTC medications
    patient_classes = {med_info['class'] for med_info in prescriptions.values()}

    for component, component_info in otc_components.items():
        # Get the class of the OTC component, if it has one
        component_class = component_info.get('class')
        if not component_class:
            continue

        for patient_class in patient_classes:
            interaction_pair = frozenset([patient_class, component_class])
            
            if interaction_pair in interaction_db:
                # An interaction was found. Identify the specific drugs involved.
                offending_prescription = [med for med, info in prescriptions.items() if info['class'] == patient_class][0]
                warning_text = interaction_db[interaction_pair]

                # Formulate the detailed counseling recommendation
                counseling_recommendation = (
                    "Based on an analysis of the medications, here is a key counseling recommendation:\n\n"
                    f"**Identified Interaction:**\n"
                    f"There is a significant drug interaction between Allison's prescription for '{offending_prescription}' "
                    f"and her use of '{otc_drug}', which contains '{component}'.\n\n"
                    f"  - '{offending_prescription}' is in a class of drugs called {patient_class}.\n"
                    f"  - '{component}' is in a class of drugs called {component_class}.\n\n"
                    "**Risk Explanation:**\n"
                    f"{warning_text}\n\n"
                    "**Pharmacist's Recommendation:**\n"
                    "The pharmacist should advise Allison to avoid taking Excedrin or other NSAIDs (like Ibuprofen) "
                    "while she is taking Fluoxetine. For her headache, a safer alternative would be a product containing "
                    "only Acetaminophen (the active ingredient in Tylenol)."
                )
                print(counseling_recommendation)
                return

    # If no interactions are found
    if not counseling_recommendation:
        print("No critical interactions were found, but general counseling on each medication's use and side effects is advised.")

if __name__ == "__main__":
    provide_counseling()