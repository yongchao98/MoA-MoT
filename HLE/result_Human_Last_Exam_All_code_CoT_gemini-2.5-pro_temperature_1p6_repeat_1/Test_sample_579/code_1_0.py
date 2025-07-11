import sys

def solve_clinical_case():
    """
    Analyzes a clinical vignette to determine the most likely diagnosis
    using a weighted scoring system based on key features.
    """
    print("Analyzing the clinical case based on provided features...")
    print("="*30)

    # --- Patient Data ---
    findings = {
        "locations": ["axillary folds", "inframammary folds", "inguinal regions"],
        "lesions": ["large bullae", "erythematous plaques", "purulent nodules"],
        "risk_factors": {"obesity_bmi_39": True, "smoking": True}
    }

    # --- Initialize Scores for each diagnosis ---
    # We will represent the final "equation" by showing the starting score and each addition.
    diagnoses = {
        "A. Malignant Intertrigo":      {"score": 0, "equation": "0"},
        "B. Allergic contact dermatitis": {"score": 0, "equation": "0"},
        "C. Hidradenitis Suppurativa":  {"score": 0, "equation": "0"},
        "D. Atopic dermatitis":         {"score": 0, "equation": "0"},
        "E. Psoriasis":                 {"score": 0, "equation": "0"}
    }

    def add_score(diagnosis_key, points, reason):
        """Helper function to update score and print reasoning."""
        diagnoses[diagnosis_key]["score"] += points
        diagnoses[diagnosis_key]["equation"] += f" + {points}"
        print(f"[{diagnosis_key}] {reason}: Adding {points} points.")

    # --- Scoring Logic ---
    print("\n1. Scoring based on Lesion Location (Intertriginous Folds):")
    # HS and Inverse Psoriasis are classic for these locations.
    add_score("C. Hidradenitis Suppurativa", 3, "Classic location in axillary, inframammary, and inguinal folds")
    add_score("E. Psoriasis", 2, "Inverse psoriasis occurs in folds")

    print("\n2. Scoring based on Lesion Type (Purulent Nodules):")
    # This is a hallmark finding for Hidradenitis Suppurativa.
    add_score("C. Hidradenitis Suppurativa", 5, "Purulent nodules are a cardinal feature")
    
    print("\n3. Scoring based on Risk Factors (Obesity and Smoking):")
    # Obesity and Smoking are very strong risk factors for HS.
    if findings["risk_factors"]["obesity_bmi_39"]:
        add_score("C. Hidradenitis Suppurativa", 2, "Obesity (BMI 39) is a major risk factor")
        add_score("E. Psoriasis", 1, "Obesity is a risk factor for inverse psoriasis")
    if findings["risk_factors"]["smoking"]:
        add_score("C. Hidradenitis Suppurativa", 2, "Smoking is a major risk factor")
        
    print("\n" + "="*30)
    print("Final Score Calculation:")
    print("="*30)

    # --- Print Final Scores and Equations ---
    highest_score = -1
    best_diagnosis_key = ""
    
    for key, value in diagnoses.items():
        final_score = value["score"]
        # The 'equation' string includes all the numbers that add up to the final score
        print(f"Diagnosis: {key}")
        print(f"Equation: {value['equation']} = {final_score}")
        print("-" * 20)
        
        if final_score > highest_score:
            highest_score = final_score
            best_diagnosis_key = key

    # --- Determine the final answer ---
    final_answer_letter = best_diagnosis_key.split('.')[0]

    print(f"\nConclusion: The diagnosis with the highest score is '{best_diagnosis_key}'.")
    print(f"The clinical picture of purulent nodules in intertriginous areas in a patient who is obese and smokes is classic for Hidradenitis Suppurativa.")
    
    # --- Final formatted output ---
    sys.stdout.flush() # Ensure all prints appear before the final answer
    print(f"<<<{final_answer_letter}>>>")

solve_clinical_case()