import sys
import io

# Redirect stdout to capture the output
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

def analyze_clinical_case():
    """
    Analyzes a clinical vignette to determine the most likely diagnosis
    by scoring potential diagnoses against patient findings.
    """
    # 1. Patient Data Extracted from Vignette
    patient_findings = {
        "age": 64,
        "bmi": 39,
        "history": {
            "locations": ["axillary folds", "inframammary folds", "inguinal regions"],
            "lesions": ["large bullae", "erythematous plaques", "purulent nodules"],
            "risk_factors": ["obesity (BMI 39)", "smoking (2-3 cigarettes/day for 15 years)"]
        }
    }

    # 2. Diagnosis Profiles
    diagnoses = {
        "A. Malignant Intertrigo": {
            "locations": ["intertriginous areas"], "lesions": ["plaques", "nodules", "ulceration"], "risk_factors": ["history of cancer"], "key_features": []
        },
        "B. Allergic contact dermatitis": {
            "locations": ["any area of contact"], "lesions": ["erythematous patches", "vesicles", "intense itching"], "risk_factors": ["allergen exposure"], "key_features": []
        },
        "C. Hidradenitis Suppurativa": {
            "locations": ["axillary folds", "inguinal regions", "inframammary folds", "perineal"],
            "lesions": ["inflammatory nodules", "abscesses", "purulent material", "sinus tracts"],
            "risk_factors": ["obesity", "smoking", "female sex", "family history"],
            "key_features": ["purulent nodules", "abscesses"]
        },
        "D. Atopic dermatitis": {
            "locations": ["flexural areas"], "lesions": ["erythematous patches", "scaling", "intense itching"], "risk_factors": ["personal/family history of atopy"], "key_features": []
        },
        "E. Psoriasis": {
            "locations": ["intertriginous areas (inverse type)"], "lesions": ["smooth erythematous plaques"], "risk_factors": ["family history", "autoimmune"], "key_features": []
        }
    }

    # 3. Scoring and Analysis
    scores = {}
    best_diagnosis = ""
    max_score = -1

    print("Analyzing patient findings against possible diagnoses...\n")

    for diagnosis, profile in diagnoses.items():
        score = 0
        reasoning = []
        
        # Check location matches
        location_matches = [loc for loc in patient_findings["history"]["locations"] if loc in profile["locations"]]
        if location_matches:
            score += len(location_matches)
            reasoning.append(f"Matches locations: {', '.join(location_matches)} ({len(location_matches)} pts)")

        # Check lesion matches
        for lesion in patient_findings["history"]["lesions"]:
            # Use 'in' for partial matches (e.g., 'purulent' in patient's "purulent nodules")
            if any(p_lesion in lesion for p_lesion in profile["lesions"]):
                score += 2  # Lesions are significant identifiers
                reasoning.append(f"Matches lesion type: '{lesion}' (2 pts)")
            # Give bonus for key features, which are highly specific
            if any(k_feature in lesion for k_feature in profile["key_features"]):
                score += 3  # Extra points for hallmark signs
                reasoning.append(f"Matches KEY feature: '{lesion}' (3 bonus pts)")

        # Check risk factor matches
        for factor in patient_findings["history"]["risk_factors"]:
            if any(p_factor in factor for p_factor in profile["risk_factors"]):
                score += 1
                reasoning.append(f"Matches risk factor: '{factor.split(' (')[0]}' (1 pt)")

        scores[diagnosis] = score
        print(f"--- Diagnosis: {diagnosis} ---")
        print(f"Final Score: {score}")
        if reasoning:
            for r in reasoning:
                print(f"  - {r}")
        else:
            print("  - Poor match based on patient's presentation.")
        print("-" * (len(diagnosis) + 20) + "\n")

        if score > max_score:
            max_score = score
            best_diagnosis = diagnosis

    print("="*50)
    print("FINAL ANALYSIS:")
    print(f"The diagnosis with the highest score is '{best_diagnosis}' with a score of {max_score}.")
    print("\nConclusion:")
    print(f"The patient is a 64-year-old woman with a BMI of {patient_findings['bmi']} (obesity) and a {patient_findings['history']['risk_factors'][1]}.")
    print("Her presentation includes lesions in three classic intertriginous areas affected by Hidradenitis Suppurativa: the axillary, inframammary, and inguinal folds.")
    print("Critically, the finding of 'purulent nodules' is a hallmark sign of Hidradenitis Suppurativa. This diagnosis best explains the combination of inflammatory lesions across multiple apocrine gland-bearing areas in a patient with significant risk factors like obesity and smoking.")
    
# Run the analysis
analyze_clinical_case()

# Restore stdout and print the captured output
sys.stdout = old_stdout
print(captured_output.getvalue())
<<<C>>>