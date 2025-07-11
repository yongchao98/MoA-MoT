def find_best_treatment():
    """
    Analyzes treatment options for a patient with fibromyalgia-like symptoms
    by scoring each option based on the symptoms it addresses.
    """
    # Key symptoms reported by the patient
    patient_symptoms = {
        "widespread_pain",
        "anxiety_and_depression",
        "sleep_issues",
        "restless_leg_syndrome",
        "paresthesia"
    }

    # Efficacy of each medication for the key symptoms
    drug_efficacy = {
        "Duloxetine": ["widespread_pain", "anxiety_and_depression"],
        "Gabapentin": ["widespread_pain", "sleep_issues", "restless_leg_syndrome", "paresthesia"],
        "cyclobenzaprine": ["sleep_issues"],
        "acetaminophen": ["widespread_pain"] # Mild efficacy, less than others
    }

    # Answer choices with the corresponding medications
    answer_choices = {
        'A': ["Duloxetine", "Gabapentin"],
        'B': ["Gabapentin"],
        'C': ["Duloxetine"],
        'D': ["cyclobenzaprine"],
        'E': ["Duloxetine", "acetaminophen"],
        'F': ["Duloxetine", "cyclobenzaprine"]
    }

    best_option = ''
    max_score = 0
    results = {}

    print("Analyzing treatment options by scoring their effectiveness against patient's symptoms...\n")

    for option, drugs in answer_choices.items():
        # Determine the set of unique symptoms covered by the drugs in the option
        symptoms_covered = set()
        for drug in drugs:
            symptoms_covered.update(drug_efficacy.get(drug, []))

        # Score is the number of unique patient symptoms covered by the treatment option
        score = len(symptoms_covered.intersection(patient_symptoms))
        results[option] = score

        if score > max_score:
            max_score = score
            best_option = option

    # Print the scoring for each option as an "equation"
    print("Scoring Results:")
    for option, score in results.items():
        drugs_in_option = " + ".join(answer_choices[option])
        # The final "equation" showing the score for each choice
        print(f"Option {option} ({drugs_in_option}) score = {score}")

    print(f"\nThe best option is '{best_option}' with a score of {max_score}, as it covers the most symptoms.")

# Run the analysis
find_best_treatment()