import collections

def analyze_clinical_case():
    """
    Analyzes a clinical case to identify the most relevant anatomical structure.
    """

    # Step 1 & 2: Analyze symptoms and map them to affected cranial nerves (CN).
    # The thoracic mass is a key localizing sign for the hoarseness.
    symptoms_to_nerves = {
        "Left facial weakness (can't lift eyebrow)": "Facial Nerve (CN VII)",
        "Loss of acoustic reflex (left ear)": "Facial Nerve (CN VII)",
        "Hoarseness & Cough (linked to thoracic mass)": "Vagus Nerve (CN X, specifically the recurrent laryngeal branch)"
    }

    print("--- Clinical Analysis ---")
    print("Patient presentation suggests involvement of the following nerves:")
    affected_nerves = set(symptoms_to_nerves.values())
    for i, (symptom, nerve) in enumerate(symptoms_to_nerves.items()):
        print(f"Symptom {i+1}: {symptom} -> Points to {nerve} lesion.")

    # Step 3: Define answer choices and their innervating nerves.
    answer_choices = {
        "A": {"structure": "Tensor tympani", "nerve": "Trigeminal Nerve (CN V)"},
        "B": {"structure": "Lateral rectus", "nerve": "Abducens Nerve (CN VI)"},
        "C": {"structure": "Intercostal muscles", "nerve": "Intercostal Nerves (Spinal)"},
        "D": {"structure": "Cricothyroid", "nerve": "Vagus Nerve (CN X)"},
        "E": {"structure": "Stylopharyngeus", "nerve": "Glossopharyngeal Nerve (CN IX)"}
    }

    print("\n--- Evaluating Answer Choices ---")
    print("Each choice is a muscle innervated by a specific nerve:")
    for key, data in answer_choices.items():
        print(f"Choice {key}: {data['structure']:<20} -> Innervated by {data['nerve']}")

    # Step 4: Find the best match between affected nerves and answer choices.
    print("\n--- Conclusion ---")
    best_match = None
    for key, data in answer_choices.items():
        if data['nerve'] in affected_nerves:
            best_match = key
            break

    if best_match:
        conclusion_structure = answer_choices[best_match]['structure']
        conclusion_nerve = answer_choices[best_match]['nerve']
        print(f"The patient's hoarseness and the thoracic mass strongly point to a lesion of the {conclusion_nerve}.")
        print(f"The {conclusion_structure} muscle (Choice {best_match}) is innervated by the {conclusion_nerve} and is directly involved in phonation (voice).")
        print("While the Facial Nerve (CN VII) is also affected, the only answer choice that matches one of the implicated cranial nerves is D.")
        print("\nTherefore, the most important anatomical structure to consider from the choices is the Cricothyroid.")

if __name__ == '__main__':
    analyze_clinical_case()