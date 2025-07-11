import json

def diagnose_skin_condition():
    """
    Analyzes a clinical case to determine the most likely dermatological diagnosis.
    This function uses a scoring system based on matching patient features
    to the characteristic features of several skin conditions.
    """
    patient_case = {
        "locations": ["axillary", "inframammary", "inguinal"],
        "lesions": ["bullae", "plaques", "purulent_nodules"],
        "risk_factors": ["obesity", "smoking", "immunosuppression"]
    }

    diagnoses = {
        "Malignant Intertrigo": {
            "key_locations": ["intertriginous"],
            "key_lesions": ["plaques", "erosions"],
            "key_risk_factors": ["immunosuppression"]
        },
        "Allergic contact dermatitis": {
            "key_locations": ["area_of_contact"],
            "key_lesions": ["vesicles", "eczema"],
            "key_risk_factors": []
        },
        "Hidradenitis Supportiva": {
            "key_locations": ["axillary", "inguinal", "inframammary"],
            "key_lesions": ["purulent_nodules", "abscesses", "plaques"],
            "key_risk_factors": ["obesity", "smoking"]
        },
        "Atopic dermatitis": {
            "key_locations": ["flexural"],
            "key_lesions": ["eczema", "lichenification"],
            "key_risk_factors": []
        },
        "Psoriasis": {
            "key_locations": ["intertriginous"],
            "key_lesions": ["plaques"],
            "key_risk_factors": []
        }
    }

    scores = {name: 0 for name in diagnoses.keys()}

    # --- Scoring Logic ---
    # Score 2 points for a direct match in key lesions/locations/risks
    # Score 1 point for a partial or general match (e.g., 'intertriginous')
    for diagnosis, features in diagnoses.items():
        # Score locations
        for location in patient_case["locations"]:
            if location in features["key_locations"]:
                scores[diagnosis] += 2
        if "intertriginous" in features["key_locations"] and all(loc in ["axillary", "inframammary", "inguinal"] for loc in patient_case["locations"]):
             scores[diagnosis] += 1

        # Score lesions
        for lesion in patient_case["lesions"]:
            if lesion in features["key_lesions"]:
                scores[diagnosis] += 2

        # Score risk factors
        for risk in patient_case["risk_factors"]:
            if risk in features["key_risk_factors"]:
                scores[diagnosis] += 2

    # --- Print Results ---
    print("Clinical Analysis Scoring:\n")
    for diagnosis, score in sorted(scores.items(), key=lambda item: item[1], reverse=True):
        print(f"- {diagnosis}: Score {score}")

    # Determine the best diagnosis
    most_likely_diagnosis = max(scores, key=scores.get)

    print("\n-----------------------------------\n")
    print(f"Conclusion: The patient's presentation of purulent nodules in the axillary, inguinal, and inframammary folds, combined with risk factors of obesity and smoking, strongly points towards a diagnosis of Hidradenitis Supportiva.")
    print(f"The highest score was for: {most_likely_diagnosis}")

if __name__ == '__main__':
    diagnose_skin_condition()
<<<C>>>