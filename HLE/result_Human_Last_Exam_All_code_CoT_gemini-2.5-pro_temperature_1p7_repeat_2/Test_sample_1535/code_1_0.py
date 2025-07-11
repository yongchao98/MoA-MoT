def solve_clinical_case():
    """
    Analyzes the clinical vignette to identify the correct anatomical region for a rash.
    This script breaks down the reasoning process based on the provided medical information.
    """

    print("Step 1: Identify Key Clinical Findings")
    key_symptoms = ["muscle weakness", "myalgia", "periorbital erythema"]
    print(f"The most significant findings are: {', '.join(key_symptoms)}.")
    
    print("\nStep 2: Formulate a Likely Diagnosis")
    print("The combination of muscle inflammation (indicated by weakness and myalgia) and skin manifestations (dermatitis) strongly suggests the diagnosis of Dermatomyositis.")
    
    print("\nStep 3: Connect Physical Exam Findings to a Specific Sign")
    exam_finding = "periorbital erythema"
    diagnosis_sign = "Heliotrope rash"
    print(f"The specific finding of '{exam_finding}' is a classic, pathognomonic sign of Dermatomyositis known as the '{diagnosis_sign}'.")
    
    print("\nStep 4: Determine the Anatomical Location of the Sign")
    anatomical_location = "Eyelids"
    print(f"The Heliotrope rash is characteristically located on the {anatomical_location}.")
    
    print("\nStep 5: Select the Corresponding Answer Choice")
    answer_choices = {
        'A': 'Dorsum of the hands',
        'B': 'Nose',
        'C': 'Eyelids',
        'D': 'Groin',
        'E': 'Shoulders'
    }
    correct_choice = 'C'
    print(f"Based on the analysis, the correct choice is the one that matches the location of the Heliotrope rash.")
    print(f"Answer Choice '{correct_choice}' corresponds to '{answer_choices[correct_choice]}'.")
    
    print("\n--- FINAL ANSWER ---")
    print("The anatomical region expected to have a rash is the Eyelids.")
    print(f"The final answer is: {correct_choice}")

# Execute the analysis function
solve_clinical_case()