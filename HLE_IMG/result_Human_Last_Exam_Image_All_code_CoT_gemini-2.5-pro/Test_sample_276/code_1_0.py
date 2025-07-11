def diagnose_patient():
    """
    This script simulates a diagnostic scoring process for a clinical case.
    It assigns points to potential diagnoses based on how well they match
    the patient's symptoms, history, and test results.
    """

    # Patient Data
    patient_features = {
        "age": "60s",
        "duration": "chronic (1 month)",
        "symptoms": ["fever", "weight loss", "diarrhea", "RLQ pain"],
        "pmh": ["uveitis", "arthritis"],
        "labs": ["leukocytosis", "fobt_positive"],
        "ct_findings": "ileocecal_thickening"
    }

    # Diagnostic criteria and scoring
    diagnoses = {
        "A. Crohn's Disease": {
            "weight": 0,
            "reasons": [],
            "rules": {
                "pmh_uveitis_arthritis": 5, # Strongest clue
                "ct_ileocecal_thickening": 3,
                "duration_chronic": 2,
                "symptoms_match": 2,
                "labs_match": 2,
            }
        },
        "C. Ileocecal Tuberculosis": {
            "weight": 0,
            "reasons": [],
            "rules": {
                "ct_ileocecal_thickening": 3, # Great mimic
                "duration_chronic": 2,
                "symptoms_match": 2,
            }
        },
        "K. Gastrointestinal Lymphoma": {
            "weight": 0,
            "reasons": [],
            "rules": {
                "age_older": 2,
                "ct_ileocecal_thickening": 2,
                "duration_chronic": 2,
                "symptoms_match": 2,
            }
        },
        "B. Yersinia Colitis": {
            "weight": 0,
            "reasons": [],
            "rules": {
                "ct_ileocecal_thickening": 2,
                "symptoms_match": 2,
                "animal_exposure_risk": 2, # Mentioned in prompt
                "duration_chronic": -1, # Usually acute
            }
        },
    }

    # Scoring logic
    print("Evaluating diagnoses based on patient data:\n")
    for name, data in diagnoses.items():
        total_score = 0
        rules = data["rules"]

        # Apply scoring rules
        if "pmh_uveitis_arthritis" in rules and all(item in patient_features["pmh"] for item in ["uveitis", "arthritis"]):
            score = rules["pmh_uveitis_arthritis"]
            total_score += score
            data["reasons"].append(f"History of uveitis/arthritis (+{score})")

        if "ct_ileocecal_thickening" in rules and patient_features["ct_findings"] == "ileocecal_thickening":
            score = rules["ct_ileocecal_thickening"]
            total_score += score
            data["reasons"].append(f"CT shows ileocecal thickening (+{score})")

        if "duration_chronic" in rules:
            score = rules["duration_chronic"]
            if score > 0 and patient_features["duration"] == "chronic (1 month)":
                 total_score += score
                 data["reasons"].append(f"Chronic duration of symptoms (+{score})")
            elif score < 0 and patient_features["duration"] == "chronic (1 month)":
                 total_score += score
                 data["reasons"].append(f"Chronic duration is atypical ({score})")

        if "symptoms_match" in rules:
            score = rules["symptoms_match"]
            total_score += score
            data["reasons"].append(f"Constitutional/GI symptoms match (+{score})")
            
        if "labs_match" in rules and all(item in patient_features["labs"] for item in ["leukocytosis", "fobt_positive"]):
            score = rules["labs_match"]
            total_score += score
            data["reasons"].append(f"Inflammatory labs match (+{score})")

        if "age_older" in rules and patient_features["age"] == "60s":
            score = rules["age_older"]
            total_score += score
            data["reasons"].append(f"Patient age is a factor (+{score})")
            
        if "animal_exposure_risk" in rules:
            score = rules["animal_exposure_risk"]
            total_score += score
            data["reasons"].append(f"Animal exposure is a risk factor (+{score})")

        data["weight"] = total_score
        
        # Print results for each diagnosis
        print(f"Diagnosis: {name}")
        print(f"  Supporting Factors: {'; '.join(data['reasons'])}")
        print(f"  Total Score: {data['weight']}\n")

    # Determine the most likely diagnosis
    most_likely_diagnosis = max(diagnoses, key=lambda k: diagnoses[k]['weight'])
    
    print("--------------------------------------------------")
    print(f"Conclusion: The most likely diagnosis is {most_likely_diagnosis}")
    print("--------------------------------------------------")


if __name__ == '__main__':
    diagnose_patient()