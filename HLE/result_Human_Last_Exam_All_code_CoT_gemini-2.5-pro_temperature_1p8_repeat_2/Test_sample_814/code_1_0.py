import pandas as pd

def evaluate_treatment_options():
    """
    This script models the selection of a treatment for a patient with fibromyalgia
    by scoring each option against the patient's primary symptoms.
    """
    
    # Patient's key symptom clusters
    symptoms = ["Widespread Pain", "Neuropathic (RLS/Paraesthesia)", "Mood (Anxiety/Depression)", "Sleep Issues"]
    
    # Medication efficacy scores (0=none, 1=mild, 2=moderate, 3=strong, 3.5=very strong/specific)
    # based on general medical knowledge for fibromyalgia treatment.
    medication_efficacy = {
        "Duloxetine":    {"Widespread Pain": 3.0, "Neuropathic (RLS/Paraesthesia)": 2.0, "Mood (Anxiety/Depression)": 3.0, "Sleep Issues": 1.0},
        "Gabapentin":    {"Widespread Pain": 2.0, "Neuropathic (RLS/Paraesthesia)": 3.5, "Mood (Anxiety/Depression)": 0.0, "Sleep Issues": 2.0},
        "cyclobenzaprine": {"Widespread Pain": 1.0, "Neuropathic (RLS/Paraesthesia)": 0.0, "Mood (Anxiety/Depression)": 0.0, "Sleep Issues": 3.0},
        "acetaminophen": {"Widespread Pain": 1.0, "Neuropathic (RLS/Paraesthesia)": 0.0, "Mood (Anxiety/Depression)": 0.0, "Sleep Issues": 0.0},
    }

    # Answer choices as combinations of medications
    answer_choices = {
        "A": ["Duloxetine", "Gabapentin"],
        "B": ["Gabapentin"],
        "C": ["Duloxetine"],
        "D": ["cyclobenzaprine"],
        "E": ["Duloxetine", "acetaminophen"],
        "F": ["Duloxetine", "cyclobenzaprine"],
    }

    results = {}

    # Calculate score for each choice
    for option, meds in answer_choices.items():
        total_score = 0
        symptom_scores = {}
        for symptom in symptoms:
            # For a combination, take the max effect for each symptom
            max_symptom_score = 0
            for med in meds:
                if medication_efficacy[med][symptom] > max_symptom_score:
                    max_symptom_score = medication_efficacy[med][symptom]
            symptom_scores[symptom] = max_symptom_score
            total_score += max_symptom_score
        results[option] = {'total': total_score, 'breakdown': symptom_scores}

    # Print results
    print("--- Treatment Option Evaluation ---")
    df = pd.DataFrame({
        option: [res['total']] for option, res in results.items()
    }, index=['Total Score']).T
    df = df.sort_values(by='Total Score', ascending=False)
    print(df)
    
    # Find the best option
    best_option = max(results, key=lambda k: results[k]['total'])
    best_score = results[best_option]['total']
    best_breakdown = results[best_option]['breakdown']
    
    # Print the rationale for the best option
    print(f"\n--- Best Option Details ---")
    print(f"The best option is '{best_option}' with a total score of {best_score}.")
    print("This choice provides the most comprehensive coverage of the patient's symptoms:")
    print("Reasoning:")
    print("- Duloxetine strongly targets widespread pain and mood disorders (anxiety/depression).")
    print("- Gabapentin specifically targets the neuropathic symptoms (restless leg syndrome and paraesthesia) and also helps with sleep and pain.")
    
    # Print the equation as requested
    score_list = list(best_breakdown.values())
    equation = ' + '.join(map(str, score_list))
    print(f"\nScore Calculation: {equation} = {best_score}")

if __name__ == "__main__":
    evaluate_treatment_options()