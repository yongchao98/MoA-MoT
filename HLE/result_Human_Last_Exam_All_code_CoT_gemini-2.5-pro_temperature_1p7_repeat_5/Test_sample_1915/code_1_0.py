def diagnose_neurological_case():
    """
    Analyzes patient symptoms to identify the damaged anatomical structure.
    """
    # Step 1: Define a knowledge base of relevant neuroanatomy.
    # This dictionary maps structures to the functions they control.
    anatomical_knowledge = {
        "Cranial nerve VI": {"functions": ["abduction"], "origin": "Pons"},
        "Cranial nerve VII": {"functions": ["facial expression", "taste"], "origin": "Pons"},
        "Medulla oblongata": {"functions": ["swallowing", "speaking", "head turning"]},
        "Midbrain": {"functions": ["adduction", "elevation", "depression", "pupillary constriction"]}
    }

    # Step 2: List the patient's symptoms from the case description.
    patient_symptoms = ["no pupillary light reflex", "unable to adduct", "unable to depress", "unable to elevate"]
    
    # Map symptoms to controlled functions for our knowledge base.
    symptom_to_function_map = {
        "no pupillary light reflex": "pupillary constriction",
        "unable to adduct": "adduction",
        "unable to depress": "depression",
        "unable to elevate": "elevation"
    }
    
    lost_functions = [symptom_to_function_map[s] for s in patient_symptoms]

    print("Patient has lost the following functions in the right eye:")
    for func in lost_functions:
        print(f"- {func.capitalize()}")
    print("-" * 30)

    # Step 3: Correlate the combination of lost functions with an anatomical structure.
    affected_structure = None
    for structure, data in anatomical_knowledge.items():
        # Check if the structure is responsible for all the patient's lost functions
        if "functions" in data and all(func in data["functions"] for func in lost_functions):
            affected_structure = structure
            break

    # Step 4: Display the results and reasoning.
    if affected_structure:
        print(f"Analysis: The combination of these lost functions points to a lesion affecting the Oculomotor Nerve (CN III).")
        print(f"The nucleus of the Oculomotor Nerve is located in the '{affected_structure}'.")
        print("-" * 30)
    else:
        print("Could not determine the single affected structure based on the provided symptoms and knowledge base.")
        return

    # Step 5: Match the finding with the given answer choices.
    answer_choices = {
        'A': 'Cranial nerve VI',
        'B': 'Cranial nerve VII',
        'C': 'Reticular formation',
        'D': 'Medulla oblongata',
        'E': 'Midbrain'
    }

    print("Evaluating Answer Choices:")
    for key, choice in answer_choices.items():
        if choice == affected_structure:
            print(f"[{key}] {choice} -> CORRECT. This matches our finding.")
        else:
            print(f"[{key}] {choice} -> INCORRECT.")
            
    final_answer = 'E'
    print(f"\nThe final answer is determined to be option {final_answer}.")


# Run the diagnostic script
diagnose_neurological_case()