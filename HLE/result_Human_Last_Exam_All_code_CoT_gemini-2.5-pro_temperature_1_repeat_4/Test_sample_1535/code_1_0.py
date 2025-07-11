def diagnose_rash_location(findings):
    """
    Analyzes clinical findings to suggest the location of a characteristic rash.
    """
    # Knowledge base mapping key symptoms to conditions and rash locations.
    conditions = {
        "Dermatomyositis": {
            "key_symptoms": {"muscle weakness", "periorbital erythema"},
            "rash_location": "Eyelids",
            "rash_name": "Heliotrope rash",
            "answer_choice": "C"
        },
        "Systemic Lupus Erythematosus": {
            "key_symptoms": {"arthralgia", "malar rash"},
            "rash_location": "Nose",
            "rash_name": "Malar (butterfly) rash",
            "answer_choice": "B"
        },
        "Psoriatic Arthritis": {
            "key_symptoms": {"arthralgia", "psoriasis"},
            "rash_location": "Dorsum of the hands",
            "rash_name": "Psoriatic plaques",
            "answer_choice": "A"
        }
    }

    # Extract key findings from the patient's case
    patient_findings = set(findings)
    best_match = None

    # Find the condition that best matches the patient's key findings
    for condition, data in conditions.items():
        if data["key_symptoms"].issubset(patient_findings):
            best_match = {
                "condition": condition,
                "location": data["rash_location"],
                "rash_name": data["rash_name"],
                "choice": data["answer_choice"]
            }
            # Dermatomyositis is the strongest match given the specific combination of symptoms
            break

    if best_match:
        print("Analysis:")
        print(f"The patient's key symptoms include muscle weakness and periorbital erythema.")
        print(f"This combination is highly suggestive of the condition: {best_match['condition']}.")
        print(f"The characteristic rash associated with this condition is the '{best_match['rash_name']}'.")
        print(f"This rash is located on the {best_match['location']}.")
        print("\nConclusion:")
        print(f"The expected anatomical region for the rash is the {best_match['location']}.")
        print(f"This corresponds to answer choice {best_match['choice']}.")
    else:
        print("The combination of symptoms does not point to a single clear diagnosis in the knowledge base.")

# Patient's key findings from the case description
patient_symptoms = ["muscle weakness", "arthralgia", "myalgia", "periorbital erythema"]

diagnose_rash_location(patient_symptoms)
