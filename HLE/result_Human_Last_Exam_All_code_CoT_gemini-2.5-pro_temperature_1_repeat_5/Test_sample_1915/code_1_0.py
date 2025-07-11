def solve_medical_case():
    """
    Analyzes patient symptoms to identify the damaged anatomical structure.
    """
    # Define the patient's observed deficits based on the clinical presentation.
    # The patient cannot adduct, depress, or elevate the right eye, and has no pupillary light reflex.
    # These functions are all associated with Cranial Nerve III.
    patient_deficits = {
        "adduction",
        "depression",
        "elevation",
        "pupillary light reflex"
    }

    # Map anatomical structures from the answer choices to the functions they control or deficits caused by their damage.
    # The midbrain is the origin of Cranial Nerve III (oculomotor nerve).
    structure_symptoms_map = {
        "A. Cranial nerve VI": {"abduction"},
        "B. Cranial nerve VII": {"facial expression"},
        "C. Reticular formation": {"consciousness"},
        "D. Medulla oblongata": {"functions of CN IX, X, XI, XII (e.g., swallowing, speech)"},
        "E. Midbrain": {"adduction", "depression", "elevation", "pupillary light reflex"}
    }

    print("Analyzing patient's symptoms:")
    print(f"Deficits: {', '.join(sorted(list(patient_deficits)))}.")
    print("\nMatching deficits to potential anatomical locations of injury:")

    best_match = None
    highest_score = 0

    # Iterate through the possible structures to find the one that explains all symptoms.
    for structure, symptoms in structure_symptoms_map.items():
        # The score is the number of patient deficits explained by damage to the structure.
        score = len(patient_deficits.intersection(symptoms))
        print(f"- Checking {structure}: Explains {score} out of {len(patient_deficits)} symptoms.")
        if score > highest_score:
            highest_score = score
            best_match = structure

    print("\n--- Conclusion ---")
    print("The patient's signs (inability to adduct, elevate, depress the eye, and loss of pupillary reflex) constitute a complete right oculomotor nerve (CN III) palsy.")
    print("The oculomotor nerve originates in the midbrain.")
    print(f"The structure that fully explains all the patient's symptoms is the {best_match.split('. ')[1]}.")
    print(f"Therefore, the correct answer is: {best_match.split('.')[0]}")


solve_medical_case()