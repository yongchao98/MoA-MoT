import sys

def solve_medical_case():
    """
    Analyzes a clinical case by scoring potential diagnoses against patient findings.
    """
    # Patient's key information from the case description
    patient_findings = {
        "age": 64,
        "bmi": 39,  # Obesity
        "habits": ["smoking"],
        "comorbidities": ["type 2 diabetes mellitus"],
        "locations": ["inframammary folds", "axillary folds", "inguinal regions"],
        "lesions": ["large bullae", "erythematous skin regions with plaques", "purulent nodules"]
    }

    # Answer choices
    diagnoses = {
        "A": "Malignant Intertrigo",
        "B": "Allergic contact dermatitis",
        "C": "Hidradenitis Suppurativa",
        "D": "Atopic dermatitis",
        "E": "Psoriasis"
    }

    # Store scores and reasoning for each diagnosis
    analysis_results = {}
    
    print("Analyzing the clinical case based on the provided findings...\n")

    # --- Evaluation Logic ---

    # A. Malignant Intertrigo
    score_A = 0
    reasoning_A = []
    reasoning_A.append("- 'Malignant Intertrigo' is not a standard diagnosis.")
    reasoning_A.append("- The presence of multiple, distinct lesion types (bullae, plaques, nodules) across multiple sites makes a single primary malignancy unlikely.")
    reasoning_A.append("- Diagnosis would require a biopsy, which is not mentioned.")
    analysis_results["A"] = {"score": score_A, "reasoning": reasoning_A}

    # B. Allergic contact dermatitis
    score_B = 0
    reasoning_B = []
    reasoning_B.append("- Location in skin folds is possible.")
    reasoning_B.append("- However, deep, purulent nodules and large bullae are not typical features of allergic contact dermatitis, which usually presents with vesicles and eczematous changes.")
    analysis_results["B"] = {"score": score_B, "reasoning": reasoning_B}

    # C. Hidradenitis Suppurativa (HS)
    score_C = 0
    reasoning_C = []
    # Check locations: HS classically affects apocrine gland-rich areas.
    if "axillary folds" in patient_findings["locations"] and "inguinal regions" in patient_findings["locations"]:
        score_C += 3
        reasoning_C.append("+ Classic locations are present (axillary, inguinal).")
    # Check lesions: HS involves inflammatory nodules and abscesses.
    if "purulent nodules" in patient_findings["lesions"]:
        score_C += 3
        reasoning_C.append("+ 'Purulent nodules' are a hallmark feature of HS.")
    if "large bullae" in patient_findings["lesions"]:
        score_C += 1
        reasoning_C.append("+ 'Large bullae' can be interpreted as abscesses, which are common in HS.")
    # Check risk factors.
    if patient_findings["bmi"] > 30 and "smoking" in patient_findings["habits"]:
        score_C += 2
        reasoning_C.append("+ Major risk factors are present (Obesity with BMI 39, Smoking).")
    analysis_results["C"] = {"score": score_C, "reasoning": reasoning_C}

    # D. Atopic dermatitis
    score_D = 0
    reasoning_D = []
    reasoning_D.append("- The primary lesion is typically an intensely pruritic (itchy) papule, leading to lichenification (skin thickening).")
    reasoning_D.append("- Deep purulent nodules are not a feature of atopic dermatitis.")
    reasoning_D.append("- While it can affect folds, the specific presentation is not typical.")
    analysis_results["D"] = {"score": score_D, "reasoning": reasoning_D}
    
    # E. Psoriasis
    score_E = 0
    reasoning_E = []
    # Inverse psoriasis affects skin folds.
    if "inframammary folds" in patient_findings["locations"] and "axillary folds" in patient_findings["locations"]:
        score_E += 1
        reasoning_E.append("+ Locations are consistent with inverse psoriasis.")
    if "erythematous skin regions with plaques" in patient_findings["lesions"]:
        score_E += 1
        reasoning_E.append("+ 'Erythematous plaques' can be seen in psoriasis.")
    reasoning_E.append("- However, 'purulent nodules' are not characteristic of psoriasis. While pustular psoriasis exists, this overall picture is a much better fit for another diagnosis.")
    analysis_results["E"] = {"score": score_E, "reasoning": reasoning_E}
    
    # --- Print Results ---
    
    highest_score = -1
    best_fit_dx = None
    
    for key, name in diagnoses.items():
        score = analysis_results[key]["score"]
        reasoning_list = analysis_results[key]["reasoning"]
        
        if score > highest_score:
            highest_score = score
            best_fit_dx = f"{key}. {name}"

        print(f"--- Diagnosis: {key}. {name} ---")
        print(f"Likelihood Score: {score}")
        for r in reasoning_list:
            print(f"  {r}")
        print("\n")

    print("--- Conclusion ---")
    print(f"The diagnosis with the highest score and most consistent clinical picture is Hidradenitis Suppurativa.")
    print("The patient presents with the classic triad of affected locations (axillary, inguinal, inframammary), hallmark lesions (purulent nodules/abscesses), and significant risk factors (obesity, smoking).")

if __name__ == "__main__":
    solve_medical_case()
<<<C>>>