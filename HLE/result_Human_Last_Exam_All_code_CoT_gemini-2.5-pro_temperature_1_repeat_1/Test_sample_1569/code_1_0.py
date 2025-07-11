def solve_medical_case():
    """
    Analyzes a clinical vignette to determine the most likely diagnosis.
    """
    patient_data = {
        "location": "Oklahoma",
        "symptoms": ["fever", "headaches", "myalgia", "disorientation"],
        "notes": "No rash mentioned. Lyme IgM can be a false positive."
    }

    # Define characteristics of potential diseases
    disease_profiles = {
        "A. Babesia microti": {
            "common_location": ["Northeast US", "Upper Midwest US"],
            "key_symptoms": ["fever", "myalgia"],
            "notes": "Neurologic symptoms are rare. Not endemic to Oklahoma."
        },
        "B. Plasmodium": {
            "common_location": ["International travel to tropics"],
            "key_symptoms": ["fever", "headaches", "myalgia"],
            "notes": "Caused by mosquito bite, not typically associated with US camping."
        },
        "C. Borrelia burgdorferi (Lyme)": {
            "common_location": ["Northeast US", "Upper Midwest US"],
            "key_symptoms": ["fever", "headaches", "myalgia"],
            "notes": "Less common in Oklahoma. Severe neurologic symptoms early on are less typical. Lab result is a distractor."
        },
        "D. Ehrlichia": {
            "common_location": ["Southeast US", "South-central US", "Oklahoma"],
            "key_symptoms": ["fever", "headaches", "myalgia", "disorientation"],
            "notes": "Transmitted by Lone Star tick, common in Oklahoma. Classic presentation matches perfectly. Rash is uncommon."
        },
        "E. Rickettsia rickettsii (RMSF)": {
            "common_location": ["Southeast US", "Oklahoma"],
            "key_symptoms": ["fever", "headaches", "myalgia", "disorientation"],
            "notes": "Also a strong candidate, but a rash (often on palms/soles) is a classic sign, and its absence is significant."
        }
    }

    best_match = None
    highest_score = -1

    print("Evaluating diagnoses based on patient data:\n")

    for disease, profile in disease_profiles.items():
        score = 0
        reasons = []

        # Score based on location
        if patient_data["location"] in profile["common_location"]:
            score += 2
            reasons.append("Location Match (+2)")
        else:
            reasons.append("Location Mismatch (+0)")

        # Score based on core symptoms
        symptom_match_count = sum(1 for symptom in patient_data["symptoms"] if symptom in profile["key_symptoms"])
        score += symptom_match_count
        reasons.append(f"Symptom Match ({symptom_match_count}/{len(patient_data['symptoms'])}) (+" + str(symptom_match_count) + ")")

        # Bonus points for specific features
        if disease == "D. Ehrlichia" and "No rash mentioned" in patient_data["notes"]:
            score += 1
            reasons.append("Absence of Rash Fits Well (+1)")
        
        if disease == "E. Rickettsia rickettsii (RMSF)" and "No rash mentioned" in patient_data["notes"]:
            score -= 1
            reasons.append("Absence of Classic Rash (-1)")


        print(f"Diagnosis: {disease}")
        equation = " + ".join(reasons)
        print(f"Calculation: {equation} = {score} total score")
        print("-" * 30)

        if score > highest_score:
            highest_score = score
            best_match = disease

    print(f"\nConclusion: The best fit based on geography and clinical presentation is {best_match}.")

solve_medical_case()