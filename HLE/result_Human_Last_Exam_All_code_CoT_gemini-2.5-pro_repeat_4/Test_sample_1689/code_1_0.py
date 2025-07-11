def solve_medical_case():
    """
    Analyzes the clinical case and evaluates the diagnostic options.
    """
    case_summary = {
        "History": "Patient started wearing new workout clothes and developed an itchy rash.",
        "Physical Exam": "Rash on the posterior axillary folds, sparing the axillary vault.",
        "Suspected Diagnosis": "Allergic contact dermatitis from textiles."
    }

    options = {
        'A': 'Skin biopsy',
        'B': 'KOH preparation',
        'C': 'Topical steroid',
        'D': 'Patch test',
        'E': 'None of the above'
    }

    reasoning = {
        'A': "A skin biopsy is invasive and not the first-line diagnostic tool for suspected allergic contact dermatitis. It's used when the diagnosis is unclear.",
        'B': "A KOH prep tests for fungus. While fungus can cause rashes in skin folds, the specific distribution and history make an allergic cause much more likely.",
        'C': "A topical steroid is a form of treatment, not a diagnostic test. It would be used after a diagnosis is made.",
        'D': "A patch test is the gold standard for diagnosing allergic contact dermatitis. It helps to identify the specific substance (like dyes or resins in clothing) that is causing the reaction. This is the most appropriate next step to confirm the suspected diagnosis.",
        'E': "This is incorrect because a definitive diagnostic test exists among the other options."
    }
    
    correct_option = 'D'

    print("Clinical Analysis:")
    for key, value in case_summary.items():
        print(f"{key}: {value}")
    
    print("\nEvaluating Diagnostic Options:")
    for option_key, option_value in options.items():
        print(f"\nOption {option_key}: {option_value}")
        print(f"Evaluation: {reasoning[option_key]}")

    print("\n" + "="*40)
    print(f"Conclusion: The best next step in diagnosis is Option {correct_option}, the '{options[correct_option]}'.")
    print("="*40)

solve_medical_case()