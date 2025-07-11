def pharmacy_counseling():
    """
    Analyzes a patient's medications to identify potential interactions and provide a counseling recommendation.
    """
    # Define the patient's prescription medications and their relevant classes/risks.
    prescriptions = {
        'Fluoxetine': {'class': 'SSRI', 'risk': 'bleeding'},
        'Atorvastatin': {'class': 'Statin', 'risk': 'muscle pain'},
        'Junel Fe': {'class': 'Oral Contraceptive', 'risk': 'blood clot'}
    }

    # Define the patient's over-the-counter medication and its key ingredients.
    otc_medication = {
        'name': 'Excedrin',
        'ingredients': {
            'Aspirin': {'class': 'NSAID', 'risk': 'bleeding'},
            'Acetaminophen': {'class': 'Analgesic', 'risk': 'liver toxicity at high dose'},
            'Caffeine': {'class': 'Stimulant', 'risk': 'insomnia/jitters'}
        }
    }

    # Check for the specific interaction between SSRIs and NSAIDs.
    recommendation = ""
    ssri_drug = None
    nsaid_ingredient = None

    for drug_name, drug_info in prescriptions.items():
        if drug_info.get('class') == 'SSRI':
            ssri_drug = drug_name
            break

    if ssri_drug:
        for ingredient_name, ingredient_info in otc_medication['ingredients'].items():
            if ingredient_info.get('class') == 'NSAID':
                nsaid_ingredient = ingredient_name
                break

    if ssri_drug and nsaid_ingredient:
        print("Interaction Found. Generating Counseling Recommendation:\n")
        
        # The final "equation" as requested, showing the interacting components.
        final_equation = f"{ssri_drug} ({prescriptions[ssri_drug]['class']}) + {nsaid_ingredient} ({otc_medication['ingredients'][nsaid_ingredient]['class']}) = Increased Bleeding Risk"
        
        recommendation = (
            f"--- Pharmacy Counseling Recommendation ---\n\n"
            f"Patient: Allison\n\n"
            f"1.  **Primary Concern:** There is a significant interaction between your prescription '{ssri_drug}' and the '{otc_medication['name']}' you took for your headache.\n\n"
            f"2.  **The Interaction:**\n"
            f"    - Your '{ssri_drug}' is an SSRI, which can increase the risk of bleeding on its own.\n"
            f"    - '{otc_medication['name']}' contains '{nsaid_ingredient}', which is an NSAID. NSAIDs also increase the risk of bleeding, particularly in the stomach.\n\n"
            f"3.  **Equation of Risk:**\n"
            f"    {final_equation}\n\n"
            f"4.  **Recommendation:**\n"
            f"    - It is advisable to **avoid** taking products containing Aspirin or other NSAIDs (like Ibuprofen) while on {ssri_drug}.\n"
            f"    - For future headaches, a safer choice would be a product containing **only Acetaminophen** (e.g., Tylenol).\n"
            f"    - Please watch for signs of stomach bleeding, such as severe stomach pain, black/tarry stools, or vomit that looks like coffee grounds. Contact your doctor if you notice any of these symptoms.\n"
        )
        print(recommendation)
    else:
        print("No major drug-drug interactions identified based on the information provided.")

if __name__ == '__main__':
    pharmacy_counseling()