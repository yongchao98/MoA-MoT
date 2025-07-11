def analyze_medical_case():
    """
    Analyzes the patient's symptoms to determine the location of the expected rash.
    """
    # Step 1: Define key patient information from the case.
    # The number 45 is the patient's age.
    patient_age = 45
    has_muscle_weakness = True
    has_periorbital_erythema = True  # This is redness around the eyes.

    print(f"Analyzing the case of the {patient_age}-year-old patient.")
    print("Key findings considered for diagnosis:")
    print(f"- Systemic symptom: Muscle weakness = {has_muscle_weakness}")
    print(f"- Physical exam finding: Periorbital erythema = {has_periorbital_erythema}\n")

    # Step 2: Correlate findings with a probable diagnosis.
    if has_muscle_weakness and has_periorbital_erythema:
        probable_diagnosis = "Dermatomyositis"
        print(f"The combination of muscle weakness and periorbital erythema strongly suggests a diagnosis of {probable_diagnosis}.\n")

        # Step 3: Identify the rash based on the specific physical finding.
        print("Dermatomyositis has several characteristic rashes.")
        print("The finding of 'periorbital erythema' (redness around the eyes) is a key feature of the 'Heliotrope rash'.")
        print("A Heliotrope rash is characteristically located on the Eyelids.\n")

        # Step 4: Evaluate the answer choices and conclude.
        answer_choices = {
            'A': 'Dorsum of the hands',
            'B': 'Nose',
            'C': 'Eyelids',
            'D': 'Groin',
            'E': 'Shoulders'
        }
        correct_choice_key = 'C'
        correct_choice_value = answer_choices[correct_choice_key]

        print(f"Conclusion: Given the specific finding of periorbital erythema, the expected anatomical region for the rash is the '{correct_choice_value}'.")
        print(f"This corresponds to answer choice: {correct_choice_key}")

analyze_medical_case()