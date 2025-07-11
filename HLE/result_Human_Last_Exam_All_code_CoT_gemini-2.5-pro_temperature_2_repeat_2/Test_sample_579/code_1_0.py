import collections

def solve_clinical_case():
    """
    Analyzes a clinical case by scoring potential diagnoses against patient features.
    """

    # Step 1: Extract Key Clinical Features from the case.
    # The patient presents with lesions in 3 classic intertriginous areas.
    # The lesions are varied: bullae (can be abscesses), plaques, and purulent nodules.
    # The patient has 2 major risk factors: obesity and smoking.
    patient_features = {
        "locations": 3,  # Axillary, inframammary, inguinal
        "lesions": 3,      # Bullae/abscesses, plaques, purulent nodules
        "risk_factors": 2  # Obesity (BMI 39), smoking
    }

    # Step 2: Define simplified profiles for each diagnosis.
    # We assign points based on how well each disease typically matches the key features.
    disease_profiles = {
        "A. Malignant Intertrigo":    {"locations": 1, "lesions": 1, "risk_factors": 1},
        "B. Allergic contact dermatitis": {"locations": 3, "lesions": 1, "risk_factors": 0},
        "C. Hidradenitis Suppurativa":  {"locations": 3, "lesions": 3, "risk_factors": 2},
        "D. Atopic dermatitis":         {"locations": 2, "lesions": 1, "risk_factors": 0},
        "E. Psoriasis":               {"locations": 3, "lesions": 1, "risk_factors": 1},
    }
    
    print("Analyzing the clinical case by scoring each potential diagnosis...\n")

    scores = {}

    # Steps 3 & 4: Calculate and display the score for each diagnosis.
    for diagnosis, profile in disease_profiles.items():
        # Match patient features against the disease profile. Max score for a feature is the patient's value.
        loc_score = min(patient_features["locations"], profile["locations"])
        les_score = min(patient_features["lesions"], profile["lesions"])
        risk_score = min(patient_features["risk_factors"], profile["risk_factors"])
        
        total_score = loc_score + les_score + risk_score
        scores[diagnosis] = total_score
        
        print(f"Diagnosis: {diagnosis}")
        print(f"Score Equation: {loc_score} (locations) + {les_score} (lesions) + {risk_score} (risk factors) = {total_score}")
        print("-" * 50)
        
    # Step 5: Identify the best match.
    best_diagnosis = max(scores, key=scores.get)

    print(f"\nThe most likely diagnosis is the one with the highest clinical match score.")
    print(f"Highest Score: {scores[best_diagnosis]}")
    print(f"Best Match: {best_diagnosis}")

solve_clinical_case()