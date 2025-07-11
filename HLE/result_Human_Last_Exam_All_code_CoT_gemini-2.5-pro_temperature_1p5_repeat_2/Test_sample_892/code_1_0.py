import pandas as pd

def solve_diagnosis():
    """
    This function analyzes a clinical case to determine the most likely diagnosis.
    It uses a weighted scoring system based on how well each diagnosis explains the patient's findings.
    """

    # Patient's clinical findings from the vignette
    findings = {
        "Dyspnea & Chronic Cough": 1,
        "Acid Reflux": 1,
        "History of COPD": 1,
        "Vertebral Mass": 10,  # Critical finding with the highest weight
        "Elevated Creatinine": 2
    }

    # Answer choices and their relevance to the findings
    diagnoses = {
        "A. Aspiration pneumonitis": {
            "Dyspnea & Chronic Cough": True, "Acid Reflux": True, "History of COPD": False, "Vertebral Mass": False, "Elevated Creatinine": False
        },
        "B. Aspiration pneumonia": {
            "Dyspnea & Chronic Cough": True, "Acid Reflux": True, "History of COPD": False, "Vertebral Mass": False, "Elevated Creatinine": False
        },
        "C. Achalasia": {
            "Dyspnea & Chronic Cough": False, "Acid Reflux": True, "History of COPD": False, "Vertebral Mass": False, "Elevated Creatinine": False
        },
        "D. Adenocarcinoma": {
            "Dyspnea & Chronic Cough": True, "Acid Reflux": True, "History of COPD": False, "Vertebral Mass": True, "Elevated Creatinine": True
        },
        "E. COPD": {
            "Dyspnea & Chronic Cough": True, "Acid Reflux": False, "History of COPD": True, "Vertebral Mass": False, "Elevated Creatinine": False
        }
    }

    print("Analyzing patient's clinical presentation...")
    print("-" * 30)
    print("Patient Findings and their Diagnostic Weight:")
    for finding, weight in findings.items():
        print(f"- {finding}: (Weight={weight})")
    print("-" * 30)

    # Calculate scores
    scores = {}
    explanation = {}
    for dx, relevant_findings in diagnoses.items():
        score = 0
        reasons = []
        for finding, is_relevant in relevant_findings.items():
            if is_relevant:
                score += findings[finding]
                reasons.append(f"+{findings[finding]} for '{finding}'")
        scores[dx] = score
        explanation[dx] = reasons

    # Determine the best diagnosis
    best_diagnosis = max(scores, key=scores.get)

    print("Evaluating each diagnosis:")
    for dx in sorted(scores.keys()):
        print(f"\nDiagnosis: {dx}")
        print(f"Explanation: {', '.join(explanation[dx])}")
        print(f"Total Score = {scores[dx]}")

    print("\n" + "=" * 30)
    print("Conclusion:")
    print("The patient has a history of COPD, which explains the chronic cough and dyspnea.")
    print("However, the new and most alarming finding is the CT scan showing a vertebral mass.")
    print("A mass on the vertebrae is highly suggestive of metastatic cancer.")
    print("Of the choices, only Adenocarcinoma (a type of lung cancer) explains the primary respiratory symptoms AND the metastatic bone mass.")
    print(f"Therefore, the most likely diagnosis is: {best_diagnosis}")
    print("=" * 30)

solve_diagnosis()
<<<D>>>