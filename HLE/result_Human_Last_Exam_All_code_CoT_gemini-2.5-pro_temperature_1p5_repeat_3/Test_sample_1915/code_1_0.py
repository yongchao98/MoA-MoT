def diagnose_neurological_issue():
    """
    Analyzes patient symptoms to identify the damaged anatomical structure.
    """
    # Step 1: Define the patient's symptoms
    patient_symptoms = {
        "no pupillary light reflex",
        "inability to adduction",
        "inability to depress",
        "inability to elevate"
    }
    
    # Step 2: Define a knowledge base mapping anatomical structures to their functions/symptoms of damage
    knowledge_base = {
        "A. Cranial nerve VI": {
            "description": "Controls outward eye movement (abduction).",
            "associated_symptoms": {"inability to abduct"}
        },
        "B. Cranial nerve VII": {
            "description": "Controls muscles of facial expression.",
            "associated_symptoms": {"facial paralysis", "loss of taste"}
        },
        "C. Reticular formation": {
            "description": "Involved in consciousness and arousal.",
            "associated_symptoms": {"altered consciousness", "coma"}
        },
        "D. Medulla oblongata": {
            "description": "Contains nuclei for CN IX, X, XI, XII (swallowing, heart rate, breathing).",
            "associated_symptoms": {"difficulty swallowing", "irregular breathing"}
        },
        "E. Midbrain": {
            "description": "Contains nuclei for Cranial Nerve III (Oculomotor Nerve).",
            "associated_symptoms": {
                "no pupillary light reflex", 
                "inability to adduction", 
                "inability to depress", 
                "inability to elevate"
            }
        }
    }
    
    print("Analyzing patient presentation...")
    print(f"Patient Symptoms: {', '.join(sorted(list(patient_symptoms)))}\n")
    
    best_match = None
    highest_match_count = -1
    
    # Step 3: Iterate through the knowledge base to find the best explanation
    for structure, data in knowledge_base.items():
        print(f"Evaluating: {structure}")
        print(f"   Function/Location: {data['description']}")
        
        # Check how many of the patient's symptoms are explained by damage to this structure
        matched_symptoms = patient_symptoms.intersection(data["associated_symptoms"])
        match_count = len(matched_symptoms)
        
        print(f"   Symptoms explained ({match_count}/{len(patient_symptoms)}): {', '.join(sorted(list(matched_symptoms))) or 'None'}")
        
        if match_count > highest_match_count:
            highest_match_count = match_count
            best_match = structure

        print("-" * 20)

    # Step 4: Print the conclusion
    print("\n---Conclusion---")
    print(f"The patient's presentation includes a complete loss of functions controlled by Cranial Nerve III.")
    print("These functions include pupillary constriction, adduction, depression, and elevation of the eye.")
    print(f"The anatomical structure that contains the nuclei for Cranial Nerve III and fully explains all the patient's symptoms is the {best_match}.")

diagnose_neurological_issue()
<<<E>>>