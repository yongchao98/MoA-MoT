import sys

def run_diagnostic_analysis():
    """
    This function analyzes a patient's clinical data to arrive at the most likely diagnosis.
    It scores each potential diagnosis based on how well it matches the patient's key findings.
    """
    
    # The patient's key findings. The normal endoscopy with abnormal imaging is the most critical clue.
    patient_findings = {
        "Heavy Smoking & Alcohol Use": {"present": True, "relevance": "Strong risk factors for Squamous Cell Carcinoma (SCC)."},
        "Imaging (Wall Thickening/Narrowing)": {"present": True, "relevance": "Indicates an infiltrative process within the esophageal wall."},
        "Endoscopy (Normal Mucosa)": {"present": True, "relevance": "Contradicts diseases that primarily affect the esophageal lining (mucosa)."}
    }

    # Potential diagnoses and their corresponding letters.
    diagnoses = {
        "A": "Streptococcal esophagitis",
        "B": "Esophageal adenocarcinoma",
        "C": "Esophageal squamous cell carcinoma",
        "D": "GERD",
        "E": "Herpes esophagitis"
    }

    # Initialize scores and a log of how they are calculated.
    scores = {dx: 0 for dx in diagnoses.values()}
    calculation_log = {dx: "0" for dx in diagnoses.values()}

    print("Analyzing patient data to determine the most likely diagnosis...\n")

    # --- Rule 1: Risk Factors ---
    print("--- Evaluating Risk Factors (Smoking & Alcohol) ---")
    score_change = 10
    diagnosis_to_update = "Esophageal squamous cell carcinoma"
    scores[diagnosis_to_update] += score_change
    calculation_log[diagnosis_to_update] += f" + {score_change}"
    print(f"Applying score for risk factors: +{score_change} points to '{diagnosis_to_update}'.\n")

    # --- Rule 2: Imaging Findings ---
    print("--- Evaluating Imaging (Wall Thickening) ---")
    score_change = 5
    # Cancers cause significant wall thickening.
    scores["Esophageal squamous cell carcinoma"] += score_change
    calculation_log["Esophageal squamous cell carcinoma"] += f" + {score_change}"
    scores["Esophageal adenocarcinoma"] += score_change
    calculation_log["Esophageal adenocarcinoma"] += f" + {score_change}"
    print(f"Applying score for imaging: +{score_change} points to cancers as they cause wall thickening.\n")

    # --- Rule 3: Endoscopy Findings ---
    print("--- Evaluating Endoscopy (Normal Mucosa) ---")
    print("This is a critical finding. It argues against diseases visible on the esophageal surface.")
    # A normal endoscopy is highly inconsistent with infectious esophagitis or typical adenocarcinoma.
    # It is the classic presentation for an infiltrative SCC.
    adjustments = {
        "Streptococcal esophagitis": -15, # Would have plaques/ulcers.
        "Esophageal adenocarcinoma": -10,  # Typically a visible mass/ulcer from the mucosa.
        "Esophageal squamous cell carcinoma": 15, # Perfectly matches an infiltrative subtype.
        "GERD": -5,                       # Severe symptoms would likely show some mucosal erythema/erosion.
        "Herpes esophagitis": -15         # Would have classic "punched-out" ulcers.
    }
    for dx, change in adjustments.items():
        op = "+" if change > 0 else "-"
        scores[dx] += change
        calculation_log[dx] += f" {op} {abs(change)}"
    print("Applying scores based on normal endoscopy...\n")
    
    # --- Final Tally ---
    print("--- Final Score Calculation ---")
    for choice_letter, diagnosis_name in diagnoses.items():
        # This fulfills the request to show each number in the final equation.
        final_equation = f"Equation for '{choice_letter}. {diagnosis_name}': {calculation_log[diagnosis_name]} = {scores[diagnosis_name]}"
        print(final_equation)

    # Determine the highest-scoring diagnosis
    most_likely_diagnosis_name = max(scores, key=scores.get)
    final_choice_letter = [k for k, v in diagnoses.items() if v == most_likely_diagnosis_name][0]

    print("\n--- Conclusion ---")
    print(f"The most likely diagnosis is '{most_likely_diagnosis_name}'.")
    print("The combination of major risk factors (smoking, alcohol) with imaging showing wall thickening despite a normal endoscopy is highly characteristic of an infiltrative Esophageal Squamous Cell Carcinoma.")

    return final_choice_letter

if __name__ == '__main__':
    # This block is for direct execution.
    # The final answer format is handled outside this script.
    # We will just print the final choice letter.
    final_answer = run_diagnostic_analysis()
    # The final step will be to wrap the answer letter in the required format.

run_diagnostic_analysis()
<<<C>>>