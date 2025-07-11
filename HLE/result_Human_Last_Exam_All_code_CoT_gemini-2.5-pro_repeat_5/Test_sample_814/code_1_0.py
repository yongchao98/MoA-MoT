import sys
import io

# Redirect stdout to capture print output for final formatting
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

def solve_clinical_case():
    """
    This function models the clinical decision-making process to find the best treatment.
    """
    # Step 1: Define patient's key symptoms
    symptoms = ['Pain', 'Anxiety/Depression', 'Sleep Issues', 'RLS/Paraesthesia']

    # Step 2: Model medication efficacy based on clinical knowledge (Scale 0-3)
    efficacy_scores = {
        'Duloxetine':     {'Pain': 3, 'Anxiety/Depression': 3, 'Sleep Issues': 1, 'RLS/Paraesthesia': 1},
        'Gabapentin':     {'Pain': 2, 'Anxiety/Depression': 1, 'Sleep Issues': 2, 'RLS/Paraesthesia': 3},
        'Cyclobenzaprine':{'Pain': 1, 'Anxiety/Depression': 0, 'Sleep Issues': 2, 'RLS/Paraesthesia': 0},
        'Acetaminophen':  {'Pain': 1, 'Anxiety/Depression': 0, 'Sleep Issues': 0, 'RLS/Paraesthesia': 0}
    }

    # Step 3: Define the provided answer choices
    treatment_options = {
        'A': 'Duloxetine+Gabapentin',
        'B': 'Gabapentin',
        'C': 'Duloxetine',
        'D': 'Cyclobenzaprine',
        'E': 'Duloxetine+Acetaminophen',
        'F': 'Duloxetine+Cyclobenzaprine'
    }

    results = {}

    # Step 4: Calculate a "Symptom Coverage Score" for each option
    for option_key, option_value in treatment_options.items():
        medications = option_value.split('+')
        total_score = 0
        
        # This will store the individual symptom scores for the final equation printout
        symptom_scores_for_option = {}

        for symptom in symptoms:
            max_symptom_score = 0
            for med in medications:
                score = efficacy_scores.get(med, {}).get(symptom, 0)
                if score > max_symptom_score:
                    max_symptom_score = score
            total_score += max_symptom_score
            symptom_scores_for_option[symptom] = max_symptom_score

        results[option_key] = {'score': total_score, 'breakdown': symptom_scores_for_option}

    # Step 5: Print the results and identify the best option
    print("Evaluating treatment options based on symptom coverage...\n")
    best_option_key = None
    max_score = -1

    for option_key, result in sorted(results.items()):
        print(f"Option {option_key} ({treatment_options[option_key]}): Score = {result['score']}")
        if result['score'] > max_score:
            max_score = result['score']
            best_option_key = option_key

    print("\n--------------------------------------------------")
    print(f"Conclusion: The best option is '{best_option_key}' with a score of {max_score}.")
    print("This option provides the most comprehensive coverage for the patient's symptoms.")
    
    # Fulfilling the requirement to print the numbers in the final equation
    best_option_breakdown = results[best_option_key]['breakdown']
    equation_parts = [f"{symptom}({score})" for symptom, score in best_option_breakdown.items()]
    print(f"Score calculation: {' + '.join(equation_parts)} = {max_score}")
    print("--------------------------------------------------")


solve_clinical_case()

# Get the captured output and print it to the actual console
output = captured_output.getvalue()
sys.stdout = old_stdout
print(output)

<<<A>>>