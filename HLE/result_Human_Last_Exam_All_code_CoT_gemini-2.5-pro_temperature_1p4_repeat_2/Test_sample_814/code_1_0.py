import sys
# Redirect print to a string to control output strictly for the final block
import io
original_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

def solve_medical_case():
    """
    This script models the treatment selection process for a complex patient case
    suggestive of Fibromyalgia. It scores each treatment option based on its
    effectiveness against the patient's primary symptoms.
    """

    # Define the patient's key symptom clusters
    symptoms = [
        "Widespread Pain",
        "Anxiety/Depression",
        "Sleep Issues",
        "Neuropathic Pain/RLS", # Paresthesia and Restless Leg Syndrome
        "Cognitive Issues"
    ]

    # Define treatment options and their estimated effectiveness scores (0=none, 1=low, 2=med, 3=high)
    # This is a simplified model for demonstration.
    treatments = {
        "A. Duloxetine+Gabapentin": {
            "Widespread Pain": 3, "Anxiety/Depression": 3, "Sleep Issues": 3,
            "Neuropathic Pain/RLS": 3, "Cognitive Issues": 1
        },
        "B. Gabapentin": {
            "Widespread Pain": 2, "Anxiety/Depression": 1, "Sleep Issues": 3,
            "Neuropathic Pain/RLS": 3, "Cognitive Issues": 1
        },
        "C. Duloxetine": {
            "Widespread Pain": 3, "Anxiety/Depression": 3, "Sleep Issues": 1,
            "Neuropathic Pain/RLS": 2, "Cognitive Issues": 1
        },
        "D. cyclobenzaprine": {
            "Widespread Pain": 1, "Anxiety/Depression": 0, "Sleep Issues": 2,
            "Neuropathic Pain/RLS": 0, "Cognitive Issues": -1 # Can worsen fog
        },
        "E. Duloxetine+acetamophen": {
            "Widespread Pain": 3, "Anxiety/Depression": 3, "Sleep Issues": 1,
            "Neuropathic Pain/RLS": 2, "Cognitive Issues": 1
        },
        "F. Duloxetine+ cyclobenzaprine": {
            "Widespread Pain": 3, "Anxiety/Depression": 3, "Sleep Issues": 2,
            "Neuropathic Pain/RLS": 2, "Cognitive Issues": 0
        }
    }

    # Calculate scores for each treatment
    scores = {}
    for option, effects in treatments.items():
        total_score = sum(effects.get(symptom, 0) for symptom in symptoms)
        scores[option] = total_score

    # Find the best option
    best_option_name = max(scores, key=scores.get)
    best_option_score = scores[best_option_name]
    best_option_details = treatments[best_option_name]

    # Print the analysis
    print("Analysis of Treatment Options for Fibromyalgia Symptoms:")
    print("-" * 55)
    print(f"Patient Symptoms: {', '.join(symptoms)}")
    print("\nScoring each option (higher is better):")
    for option, score in sorted(scores.items(), key=lambda item: item[1], reverse=True):
        print(f"  {option:<28} | Total Score: {score}")

    print("\n--- Conclusion ---")
    print(f"The best option is '{best_option_name}' because it provides the most comprehensive coverage for the patient's specific set of symptoms.")

    # Create the final equation string as requested
    equation_numbers = [str(best_option_details.get(s, 0)) for s in symptoms]
    equation_str = " + ".join(equation_numbers)

    print("\nThe score for the best option is calculated by summing its effectiveness for each symptom:")
    print(f"Effectiveness scores: {', '.join(f'{s}: {best_option_details.get(s, 0)}' for s in symptoms)}")
    print(f"\nFinal Equation: {equation_str} = {best_option_score}")


# Execute the function and capture the output
solve_medical_case()
# Get the content from the captured output
final_output = captured_output.getvalue()
# Restore original stdout
sys.stdout = original_stdout
# Print the final captured output in the desired format
print(final_output)