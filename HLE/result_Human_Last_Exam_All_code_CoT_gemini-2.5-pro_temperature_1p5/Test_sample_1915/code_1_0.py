def solve_medical_case():
    """
    Analyzes a clinical case to identify the damaged anatomical structure.
    This script programmatically walks through the reasoning.
    """

    # --- Step 1: Define Patient's Key Information and Symptoms ---
    patient_age = 33
    right_eye_symptoms = {
        "Pupillary Light Reflex": "Absent",
        "Adduction (inward movement)": "Unable",
        "Depression (downward movement)": "Unable",
        "Elevation (upward movement)": "Unable"
    }

    print("Analyzing the case of a {}-year-old female...".format(patient_age))
    print("\nKey Symptoms in the Right Eye:")
    for symptom, status in right_eye_symptoms.items():
        print(f"- {symptom}: {status}")

    # --- Step 2: Define Functions of Relevant Anatomical Structures ---
    # This represents our medical knowledge base.
    cranial_nerve_III_functions = [
        "Adduction of the eye",
        "Elevation of the eye",
        "Depression of the eye",
        "Pupillary constriction (light reflex)"
    ]

    anatomical_structures = {
        "A. Cranial nerve VI": "Controls abduction (outward eye movement).",
        "B. Cranial nerve VII": "Controls muscles of facial expression.",
        "C. Reticular formation": "Involved in consciousness and arousal.",
        "D. Medulla oblongata": "Contains nuclei for cranial nerves IX, X, XI, XII.",
        "E. Midbrain": "Contains the nucleus for Cranial Nerve III (Oculomotor nerve)."
    }

    # --- Step 3: Perform Logical Deduction ---
    print("\n--- Logical Deduction ---")
    print("1. The patient's symptoms (impaired adduction, elevation, depression, and pupillary reflex) perfectly match a complete palsy of Cranial Nerve III.")

    print("\n2. The logical question is: Where does Cranial Nerve III originate?")

    print("\n3. Reviewing the answer choices:")
    correct_location = None
    correct_choice = None
    for choice, function in anatomical_structures.items():
        # The key to the logical 'equation' is finding which structure contains CN III.
        if "Cranial Nerve III" in function:
            correct_location = choice
            correct_choice = choice[0]
            break

    print(f"   - The knowledge base indicates that the '{correct_location}' is the location of the nucleus for Cranial Nerve III.")

    print("\n--- Final Conclusion ---")
    print("The patient's presentation is explained by a lesion affecting Cranial Nerve III.")
    # Here is the 'final equation' as a logical statement
    print("The final 'equation' is: Patient's Symptoms = Cranial Nerve III Palsy. The location of the CN III nucleus = Midbrain.")
    print(f"Therefore, damage to the '{correct_location}' explains the patient's presentation.")
    print(f"The correct answer choice is {correct_choice}.")

    # --- Final Answer in Specified Format ---
    return f"<<<{correct_choice}>>>"

# Execute the function and print the final answer
final_answer = solve_medical_case()
print(final_answer)
