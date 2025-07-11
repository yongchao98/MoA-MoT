import operator

def evaluate_fibromyalgia_treatment():
    """
    This script evaluates treatment options for a patient with symptoms of Fibromyalgia.
    It scores each option based on its effectiveness against a set of key symptoms.
    """
    # 1. Define patient's key symptoms
    symptoms = ['widespread_pain', 'fatigue', 'anxiety_depression', 'sleep_issues', 'rls_paresthesia']

    # 2. Model the efficacy of each drug for each symptom (0=none, 1=moderate, 2=high)
    # This is a simplified model for demonstration.
    drug_efficacy = {
        'Duloxetine':     {'widespread_pain': 2, 'fatigue': 1, 'anxiety_depression': 2, 'sleep_issues': 1, 'rls_paresthesia': 1},
        'Gabapentin':     {'widespread_pain': 2, 'fatigue': 0, 'anxiety_depression': 1, 'sleep_issues': 2, 'rls_paresthesia': 2},
        'cyclobenzaprine':{'widespread_pain': 1, 'fatigue': 0, 'anxiety_depression': 0, 'sleep_issues': 2, 'rls_paresthesia': 0},
        'acetaminophen':  {'widespread_pain': 0, 'fatigue': 0, 'anxiety_depression': 0, 'sleep_issues': 0, 'rls_paresthesia': 0}
    }

    # 3. Define the answer choices and the drugs they contain
    treatment_options = {
        "A. Duloxetine+Gabapentin": ["Duloxetine", "Gabapentin"],
        "B. Gabapentin": ["Gabapentin"],
        "C. Duloxetine": ["Duloxetine"],
        "D. cyclobenzaprine": ["cyclobenzaprine"],
        "E. Duloxetine+acetaminophen": ["Duloxetine", "acetaminophen"],
        "F. Duloxetine+cyclobenzaprine": ["Duloxetine", "cyclobenzaprine"]
    }

    print("Evaluating treatment options based on symptom coverage...\n")

    results = {}

    # 4. Loop through each option, calculate, and print the score
    for option, drugs in treatment_options.items():
        total_score = 0
        symptom_scores = {}
        for symptom in symptoms:
            # For combination therapy, take the score of the most effective drug for that symptom
            best_score_for_symptom = 0
            for drug in drugs:
                score = drug_efficacy[drug][symptom]
                if score > best_score_for_symptom:
                    best_score_for_symptom = score
            symptom_scores[symptom] = best_score_for_symptom
            total_score += best_score_for_symptom
        
        results[option] = {'total': total_score, 'breakdown': symptom_scores}

        print(f"Option: {option}")
        print(f"Symptom Score Breakdown: {symptom_scores}")
        print(f"Total Score: {total_score}\n")

    # 5. Determine the best option
    best_option_name = max(results, key=lambda k: results[k]['total'])
    best_option_data = results[best_option_name]

    print("---" * 10)
    print(f"Conclusion: The best option is '{best_option_name}'")
    print("It provides the most comprehensive coverage for the patient's symptoms.")
    
    # Fulfills the requirement to "output each number in the final equation"
    breakdown = best_option_data['breakdown']
    equation_parts = [f"{symptom.replace('_',' ')} ({score})" for symptom, score in breakdown.items()]
    equation = " + ".join(equation_parts)
    print(f"Final Score Calculation: {equation} = {best_option_data['total']}")


if __name__ == '__main__':
    evaluate_fibromyalgia_treatment()
