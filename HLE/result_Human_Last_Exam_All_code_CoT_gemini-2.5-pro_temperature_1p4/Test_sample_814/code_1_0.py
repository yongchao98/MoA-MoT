def solve_medical_question():
    """
    This function models the decision-making process for selecting the best treatment
    based on symptom coverage.
    """

    # Define the patient's key symptoms that need treatment
    patient_symptoms = {
        'widespread_pain',
        'anxiety_depression',
        'sleep_issues',
        'paresthesia',  # Neuropathic pain component
        'restless_leg_syndrome'
    }

    # Define what each medication primarily treats, based on clinical knowledge
    treatment_coverage = {
        'Duloxetine': {'widespread_pain', 'anxiety_depression', 'paresthesia'},
        'Gabapentin': {'widespread_pain', 'paresthesia', 'sleep_issues', 'restless_leg_syndrome'},
        'Cyclobenzaprine': {'sleep_issues'}, # Primarily for sleep at low doses for FM
        'Acetaminophen': set() # Not a primary treatment for any of these key FM symptoms
    }

    # Define the answer choices as combinations of the medications
    answer_choices = {
        'A': ['Duloxetine', 'Gabapentin'],
        'B': ['Gabapentin'],
        'C': ['Duloxetine'],
        'D': ['Cyclobenzaprine'],
        'E': ['Duloxetine', 'Acetaminophen'],
        'F': ['Duloxetine', 'Cyclobenzaprine']
    }

    # Calculate a "suitability score" for each choice
    best_choice = ''
    highest_score = -1

    print("Evaluating treatment options based on symptom coverage:\n")
    print(f"Patient's Key Symptoms: {', '.join(sorted(list(patient_symptoms)))}\n")

    for option, medications in answer_choices.items():
        # Combine the coverage of all meds in the option
        combined_coverage = set()
        for med in medications:
            combined_coverage.update(treatment_coverage[med])
        
        # Calculate score as the number of patient symptoms covered
        covered_symptoms = patient_symptoms.intersection(combined_coverage)
        score = len(covered_symptoms)

        print(f"Option {option} ({' + '.join(medications)}):")
        print(f"  - Covers {score} symptoms: {', '.join(sorted(list(covered_symptoms)))}")
        print(f"  - Score = {score}\n")
        
        if score > highest_score:
            highest_score = score
            best_choice = option

    print("-----------------------------------------")
    print(f"Conclusion: Option {best_choice} provides the most comprehensive coverage of the patient's symptoms.")

solve_medical_question()