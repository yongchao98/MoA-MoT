import textwrap

def generate_pharmacist_counseling():
    """
    Analyzes a patient's medications to identify interactions and generate counseling recommendations.
    """
    # 1. Represent patient's medications
    # The numbers (20, 1.5, 30) from her prescriptions are noted here.
    patient_meds = {
        'Fluoxetine': '20mg',
        'Junel Fe': '1.5/30mg',
        'Atorvastatin': '20mg',
        'Excedrin': 'Unknown (OTC)'
    }

    # 2. Drug knowledge base
    # Mapping medications to their relevant drug class for interaction checking
    drug_classes = {
        'Fluoxetine': 'SSRI',
        'Aspirin': 'NSAID',
        'Junel Fe': 'Oral Contraceptive'
    }

    # OTC medications and their relevant active ingredients
    otc_components = {
        'Excedrin': ['Aspirin', 'Acetaminophen', 'Caffeine']
    }

    # 3. Interaction knowledge base
    # Key: a tuple of interacting drug classes
    # Value: the description of the interaction
    interaction_db = {
        ('SSRI', 'NSAID'): "significantly increases the risk of gastrointestinal (GI) bleeding"
    }

    # --- Analysis Logic ---
    active_drug_classes = {}
    for med, dosage in patient_meds.items():
        # Check OTC components first
        if med in otc_components:
            for component in otc_components[med]:
                if component in drug_classes:
                    active_drug_classes[drug_classes[component]] = f"{med} (contains {component})"
        # Check prescription meds
        if med in drug_classes:
            active_drug_classes[drug_classes[med]] = f"{med} {dosage}"

    # Identify and report interactions
    recommendation_found = False
    print("Pharmacist Counseling Recommendation for Allison:")
    print("-" * 50)

    # Check for the primary interaction
    if 'SSRI' in active_drug_classes and 'NSAID' in active_drug_classes:
        ssri_med_name = active_drug_classes['SSRI']
        nsaid_med_name = active_drug_classes['NSAID']
        interaction_desc = interaction_db[('SSRI', 'NSAID')]

        print("Primary concern: A significant drug interaction was identified.")
        print(f"  - Drug 1: {ssri_med_name} (Class: SSRI)")
        print(f"  - Drug 2: {nsaid_med_name} (Class: NSAID)")
        print(f"\nInteraction: Taking these two medications together {interaction_desc}.")

        print("\nRecommendation:")
        rec_text = (f"1. You should avoid taking Excedrin or other products containing Aspirin or NSAIDs "
                    f"(like ibuprofen) for your headache while on {ssri_med_name}. "
                    "For headache relief, a safer alternative would be a product with only acetaminophen (like Tylenol). "
                    "Please discuss this with your doctor or me before taking any new medicine.\n\n"
                    "2. Please be aware of the signs of GI bleeding, such as unusual bruising, "
                    "black or tarry stools, or vomit that looks like coffee grounds. "
                    "If you notice any of these, contact your doctor immediately.")
        print(textwrap.fill(rec_text, width=70))
        recommendation_found = True

    # Add secondary counseling point about headache and oral contraceptives
    if 'Oral Contraceptive' in active_drug_classes:
        oc_med_name = active_drug_classes['Oral Contraceptive']
        print("\n" + "-"*50)
        print("Secondary Point:")
        rec2_text = (f"Additionally, since you are experiencing headaches while taking {oc_med_name}, "
                     "it is a good idea to mention this to your doctor. They will want to ensure it is not related "
                     "to your birth control medication.")
        print(textwrap.fill(rec2_text, width=70))


    if not recommendation_found:
        print("No critical interactions found with the provided information, but it's always best to review all medications with your pharmacist.")

if __name__ == "__main__":
    generate_pharmacist_counseling()