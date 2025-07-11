def solve_medical_case():
    """
    This script analyzes the patient's symptoms to find the correct answer.
    """
    # Step 1: Define the key signs and symptoms from the case description.
    key_symptoms = ["muscle weakness", "myalgia", "arthralgia"]
    key_physical_signs = ["periorbital erythema"]

    print(f"Patient's key findings: {', '.join(key_symptoms)}, and {key_physical_signs[0]}.")

    # Step 2: Correlate findings to a likely diagnosis.
    diagnosis = "Dermatomyositis"
    print(f"The combination of muscle weakness and periorbital erythema (a 'heliotrope rash') strongly suggests a diagnosis of {diagnosis}.")

    # Step 3: Identify other classic signs of the diagnosis.
    # Dermatomyositis has two hallmark skin findings.
    heliotrope_rash = "A rash on the Eyelids (periorbital region)."
    gottrons_sign = "A rash on the Dorsum of the hands (over the knuckles)."
    print(f"Dermatomyositis is characterized by two main rashes: {heliotrope_rash} and {gottrons_sign}.")

    # Step 4: Determine the answer.
    print(f"The patient already has periorbital erythema, which is the {heliotrope_rash}.")
    print(f"Therefore, the other expected rash is Gottron's sign, which appears on the Dorsum of the hands.")

    # Step 5: Match with the given answer choices.
    answer_choices = {
        "A": "Dorsum of the hands",
        "B": "Nose",
        "C": "Eyelids",
        "D": "Groin",
        "E": "Shoulders"
    }
    correct_answer_letter = "A"
    print("\nMatching this location to the answer choices:")
    print(f"Choice {correct_answer_letter} is '{answer_choices[correct_answer_letter]}'.")

solve_medical_case()