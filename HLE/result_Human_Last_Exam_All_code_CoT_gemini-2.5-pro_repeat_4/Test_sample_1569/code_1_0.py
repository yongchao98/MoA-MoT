def solve_medical_case():
    """
    This function analyzes a clinical vignette to determine the most likely diagnosis
    by scoring potential diseases based on key patient data.
    """
    # Patient Data from the vignette
    patient_data = {
        "location": "Oklahoma",
        "symptoms": ["fever", "headaches", "myalgia", "disorientation", "heart_murmur"],
        "labs": "Positive Lyme IgM"
    }

    # Database of disease characteristics
    diseases = {
        "A. Babesia microti": {
            "geography": ["Northeast US", "Upper Midwest US"],
            "core_symptoms": ["fever", "headaches", "myalgia"],
            "severe_symptoms": []
        },
        "B. Plasmodium": {
            "geography": ["Malaria-endemic regions (not US)"],
            "core_symptoms": ["fever", "headaches", "myalgia"],
            "severe_symptoms": ["disorientation"]
        },
        "C. Borrelia burgdorferi": {
            "geography": ["Northeast US", "Upper Midwest US"],
            "core_symptoms": ["fever", "headaches", "myalgia"],
            "severe_symptoms": ["disorientation", "heart_murmur"]
        },
        "D. Ehrlichia": {
            "geography": ["South-central US", "Southeast US", "Oklahoma"],
            "core_symptoms": ["fever", "headaches", "myalgia"],
            "severe_symptoms": ["disorientation"]
        },
        "E. Rickettsia rickettsii": {
            "geography": ["South-central US", "Southeast US", "Oklahoma"],
            "core_symptoms": ["fever", "headaches", "myalgia"],
            "severe_symptoms": ["disorientation", "heart_murmur"]
        }
    }

    scores = {}
    print("Scoring each potential diagnosis based on patient data:\n")

    for disease, data in diseases.items():
        score = 0
        
        # Scoring based on geography (crucial factor)
        # We give a high weight to a geographic match.
        geo_score = 3 if patient_data["location"] in data["geography"] else 0
        score += geo_score

        # Scoring based on symptom match
        symptom_score = 0
        for symptom in patient_data["symptoms"]:
            if symptom in data["core_symptoms"] or symptom in data["severe_symptoms"]:
                symptom_score += 1
        score += symptom_score

        scores[disease] = score
        
        # The 'equation' is the sum of scores for different clinical features.
        # We print each component to show the calculation.
        print(f"Diagnosis: {disease}")
        print(f"Calculation: Geographic Match Score ({geo_score}) + Symptom Match Score ({symptom_score})")
        print(f"Final Score = {score}\n")
        
    # Find the disease with the highest score
    most_likely_diagnosis = max(scores, key=scores.get)
    
    print("---Conclusion---")
    print(f"The most likely diagnosis is the one with the highest score: {most_likely_diagnosis}")

solve_medical_case()
<<<E>>>