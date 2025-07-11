import collections

def diagnose_neurological_lesion():
    """
    Analyzes patient symptoms to identify the location of a neurological lesion.
    """
    # Define the patient's key symptoms and the cranial nerve responsible
    # The patient has a complete right oculomotor nerve (CN III) palsy.
    patient_symptoms = {
        "no pupillary light reflex (motor part)": "CN III",
        "inability to adduct eye": "CN III",
        "inability to depress eye": "CN III",
        "inability to elevate eye": "CN III",
    }

    # Define the answer choices and the key neurological structures/functions they are associated with.
    # The Midbrain is the key structure as it contains the nucleus for CN III.
    answer_choices = {
        "A": {"name": "Cranial nerve VI", "functions": ["eye abduction"]},
        "B": {"name": "Cranial nerve VII", "functions": ["facial expression", "taste"]},
        "C": {"name": "Reticular formation", "functions": ["arousal", "consciousness"]},
        "D": {"name": "Medulla oblongata", "functions": ["nuclei of CN IX, X, XI, XII"]},
        "E": {"name": "Midbrain", "functions": ["nuclei of CN III, IV"]},
    }

    print("Analyzing patient presentation...")
    print("-" * 40)
    print("Symptoms Reported:")
    for symptom in patient_symptoms:
        print(f"- {symptom}")

    print("\nAnalysis:")
    print("All of the patient's ocular symptoms (impaired adduction, depression, elevation, and pupillary reflex) are caused by a deficit in the Oculomotor Nerve (Cranial Nerve III).")
    print("The task is to find the anatomical structure from the choices that is responsible for CN III function.")

    print("\nEvaluating Choices:")
    best_choice = ""
    for choice, details in answer_choices.items():
        # Check if the structure contains the nucleus for CN III
        if "nuclei of CN III, IV" in details["functions"]:
            best_choice = choice
            print(f"({choice}) {details['name']}: Contains the nucleus for Cranial Nerve III. Damage here directly explains a CN III palsy. This is a strong candidate.")
        else:
            print(f"({choice}) {details['name']}: Does not contain the CN III nucleus and cannot explain the combination of symptoms.")

    print("-" * 40)
    print("Conclusion:")
    if best_choice:
        conclusion_text = f"The damage is to the {answer_choices[best_choice]['name']}, as it contains the nucleus for the Oculomotor Nerve (CN III), explaining all the patient's eye movement and pupillary deficits."
        print(conclusion_text)
        print(f"Therefore, the correct answer is {best_choice}.")
    else:
        print("Could not determine the correct answer based on the provided data.")

# Run the diagnosis
diagnose_neurological_lesion()