def solve_medical_mystery():
    """
    Analyzes a patient's clinical case using a simple scoring model
    to determine the most likely disease.
    """
    # Define potential diseases and their key associated features.
    diseases = {
        "Granulomatosis with Polyangiitis (GPA)": [
            "older male", "smoker", "fatigue", "arthritis", "shortness of breath",
            "pulmonary nodules", "confusion", "bruising", "difficulty swallowing",
            "started on steroids", "immunosuppression", "septic shock"
        ],
        "Lung Cancer with Paraneoplastic Syndrome": [
            "older male", "smoker", "ship building", "fatigue", "loss of appetite",
            "shortness of breath", "pulmonary nodules", "difficulty swallowing",
            "confusion", "arthritis"
        ],
        "Sarcoidosis": [
            "fatigue", "arthritis", "shortness of breath", "pulmonary nodules",
            "cutaneous lesions", "started on steroids"
        ],
        "Systemic Lupus Erythematosus (SLE)": [
             "fatigue", "swelling", "pain in his wrists, ankles, and elbows", "bruising", "confusion", "cutaneous lesions", "fever"
        ],
        "Tuberculosis (TB)": [
            "fatigue", "loss of appetite", "shortness of breath", "pulmonary nodules",
            "fever", "productive cough", "immunosuppression"
        ]
    }

    # List the key data points from the patient's case.
    patient_profile = [
        "62-year-old man", "older male",
        "20-pack-year history of smoking", "smoker",
        "work in ship building", "ship building",
        "fatigue",
        "swelling, and pain in his wrists, ankles, and elbows", "arthritis",
        "confusion",
        "bruising",
        "difficulty swallowing",
        "loss of appetite",
        "shortness of breath",
        "multiple pulmonary nodules", "pulmonary nodules",
        "started on steroids",
        "fever",
        "productive cough with green sputum", "productive cough",
        "cutaneous lesions",
        "ineffective Aminoglycoside therapy",
        "died from septic shock", "septic shock"
    ]

    # Initialize a dictionary to hold the score for each disease.
    scores = {disease: 0 for disease in diseases}
    
    # --- Scoring Calculation ---
    print("--- Differential Diagnosis Analysis ---")
    print("Matching patient profile against disease features:\n")

    # Iterate through each disease and its features to calculate a score.
    for disease, features in diseases.items():
        print(f"Evaluating: {disease}")
        disease_score = 0
        matching_features = []
        for feature in features:
            # Check if any part of the patient's profile contains the feature.
            if any(feature in profile_item for profile_item in patient_profile):
                disease_score += 1
                matching_features.append(feature)
        
        scores[disease] = disease_score
        print(f"  Score: {disease_score}")
        print(f"  Matching Features: {', '.join(matching_features)}\n")

    # Determine the most likely diagnosis based on the highest score.
    most_likely_disease = max(scores, key=scores.get)

    # --- Final Conclusion ---
    print("--- Conclusion ---")
    print("\nFinal Scores:")
    # Sort diseases by score in descending order for a clear ranking.
    sorted_scores = sorted(scores.items(), key=lambda item: item[1], reverse=True)
    for disease, score in sorted_scores:
        print(f"{disease}: {score}")

    print(f"\nThe patient's history and symptoms most closely match the profile of: {most_likely_disease}")

solve_medical_mystery()