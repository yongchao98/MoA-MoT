import sys
# Redirecting stderr to null to suppress any potential non-stdout messages in the execution environment
# for a cleaner output, though none are expected here.
sys.stderr = open('/dev/null', 'w')

def analyze_patient_case():
    """
    Analyzes a clinical case using a simplified scoring model to determine the most likely diagnosis.
    This function is for illustrative purposes only and is not a substitute for professional medical advice.
    """

    # --- Patient Data ---
    # Data extracted from the case description
    patient_data = {
        "age": 57,
        "symptoms": ["dyspnea", "chronic cough", "acid reflux"],
        "history": ["COPD"],
        "findings": {
            "ct_scan": "mass of the vertebrae",
            "serum_creatinine": 2.1  # Assuming "blood urine creatine" is a typo for serum creatinine
        }
    }

    # --- Scoring Model ---
    # Assigns points to findings that support a diagnosis of metastatic adenocarcinoma
    scoring_points = {
        "COPD_risk_factor": 2,  # COPD is a major risk factor for lung cancer
        "respiratory_symptoms": 1,  # For dyspnea and cough
        "metastasis_sign": 5    # Vertebral mass is a strong indicator of metastasis
    }

    # --- Analysis ---
    # Calculate the score for Adenocarcinoma, which is the most likely diagnosis
    # to explain the full clinical picture.
    print("--- Diagnostic Reasoning Process ---")
    print("Patient has respiratory symptoms (dyspnea, cough) and a history of COPD.")
    print("While this could be a COPD exacerbation, a new critical finding must be explained.")
    print("Critical Finding: A mass on the vertebrae revealed by CT scan.\n")
    print("This finding is highly concerning for cancer that has spread to the bone (metastasis).")
    print("Given the patient's COPD history (a major risk factor), lung cancer is a primary suspect.\n")
    
    # Building the 'equation' based on the scoring model
    score = 0
    equation_components = []

    # Score for risk factor
    risk_factor_score = scoring_points["COPD_risk_factor"]
    score += risk_factor_score
    equation_components.append(str(risk_factor_score))
    print(f"Points for COPD as a risk factor: {risk_factor_score}")

    # Score for symptoms
    symptom_score = scoring_points["respiratory_symptoms"] * 2 # dyspnea + cough
    score += symptom_score
    equation_components.append(str(scoring_points["respiratory_symptoms"]))
    equation_components.append(str(scoring_points["respiratory_symptoms"]))
    print(f"Points for respiratory symptoms (dyspnea & cough): {symptom_score}")

    # Score for key finding
    metastasis_score = scoring_points["metastasis_sign"]
    score += metastasis_score
    equation_components.append(str(metastasis_score))
    print(f"Points for vertebral mass (sign of metastasis): {metastasis_score}\n")

    # --- Final Conclusion & Equation ---
    # The other diagnoses (Aspiration, Achalasia, COPD exacerbation) fail to
    # account for the vertebral mass, making them far less likely.
    
    # We output each number in the final equation.
    print("Scoring equation for the most likely diagnosis (Adenocarcinoma):")
    final_equation = f"{equation_components[0]} + {equation_components[1]} + {equation_components[2]} + {equation_components[3]} = {score}"
    print(final_equation)

    print("\nConclusion: The combination of risk factors, symptoms, and especially the evidence of a vertebral mass makes Adenocarcinoma (likely metastatic lung cancer) the most unifying diagnosis.")


if __name__ == '__main__':
    analyze_patient_case()