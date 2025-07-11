def find_rash_location():
    """
    Analyzes a clinical case to determine the expected location of a rash
    based on pathognomonic signs.
    """
    # Step 1: Define the key clinical findings from the case.
    patient_findings = {
        "systemic": ["muscle weakness", "myalgia", "arthralgia", "fatigue"],
        "physical_exam": "periorbital erythema"
    }

    print("Analyzing the patient's clinical presentation...")
    print(f"Systemic symptoms: {', '.join(patient_findings['systemic'])}")
    print(f"Key physical exam finding: {patient_findings['physical_exam']}")
    print("-" * 30)

    # Step 2: Correlate the key finding with a specific medical sign.
    key_sign = patient_findings["physical_exam"]
    medical_term = "heliotrope rash"
    diagnosis = "Dermatomyositis"

    print(f"The finding '{key_sign}' combined with systemic muscle symptoms is highly suggestive of {diagnosis}.")
    print(f"'{key_sign}' is the clinical description for a '{medical_term}'.")
    print("-" * 30)

    # Step 3: Identify the anatomical location of this specific sign.
    rash_location = "Eyelids"
    print(f"A {medical_term} is a purplish rash located on the '{rash_location}'.")
    print("-" * 30)

    # Step 4: Compare with the given answer choices.
    answer_choices = {
        "A": "Dorsum of the hands",
        "B": "Nose",
        "C": "Eyelids",
        "D": "Groin",
        "E": "Shoulders"
    }

    print("Matching the location to the answer choices:")
    for choice, description in answer_choices.items():
        print(f"Choice {choice}: {description}")
        if description == rash_location:
            correct_choice = choice

    print("-" * 30)
    print(f"The expected anatomical region for the rash is the '{rash_location}', which corresponds to answer choice {correct_choice}.")

# Execute the diagnostic logic.
find_rash_location()