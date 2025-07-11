import collections

def solve_clinical_case():
    """
    Analyzes a clinical vignette to determine the most likely diagnosis.
    """
    # Step 1: Define patient's clinical data
    patient_symptoms = {"fever", "headaches", "myalgia", "disorientation", "heart_murmur"}
    patient_history = {"location": "Oklahoma"}
    patient_labs = {"Lyme_IgM": "positive"}

    # Step 2: Define disease profiles
    # Data represents classic presentations.
    disease_profiles = {
        "A. Babesia microti": {
            "symptoms": {"fever", "headaches", "myalgia"},
            "geography": ["Northeast US", "Upper Midwest US"],
            "notes": "Cardiac and CNS signs are uncommon."
        },
        "B. Plasmodium": {
            "symptoms": {"fever", "headaches", "myalgia"},
            "geography": [], # Not endemic in Oklahoma
            "notes": "Acquired via travel to endemic zones, not camping in Oklahoma."
        },
        "C. Borrelia burgdorferi": {
            "symptoms": {"fever", "headaches", "myalgia", "disorientation", "heart_murmur"},
            "geography": ["Northeast US", "Upper Midwest US", "Oklahoma"],
            "notes": "Classic for early disseminated disease (neuro/cardiac)."
        },
        "D. Ehrlichia": {
            "symptoms": {"fever", "headaches", "myalgia", "disorientation"},
            "geography": ["Oklahoma", "Southeastern US"],
            "notes": "High incidence in Oklahoma, but cardiac signs are rare."
        },
        "E. Rickettsia rickettsii": {
            "symptoms": {"fever", "headaches", "myalgia", "disorientation"},
            "geography": ["Oklahoma", "Southeastern US"],
            "notes": "High incidence in Oklahoma, but characteristic rash is absent."
        }
    }

    # Step 3: Develop scoring logic and calculate scores
    scores = collections.defaultdict(int)
    
    print("Diagnostic Scoring Breakdown:")
    print("="*30)
    for disease, profile in disease_profiles.items():
        # Base score for non-specific symptoms
        base_symptoms = {"fever", "headaches", "myalgia"}
        symptom_match = base_symptoms.issubset(profile["symptoms"])
        if symptom_match:
            scores[disease] += 1
            print(f"'{disease}': Base symptoms match. Score +1. Current score: {scores[disease]}")

        # Score for specific neurologic sign (disorientation)
        if "disorientation" in patient_symptoms and "disorientation" in profile["symptoms"]:
            scores[disease] += 2
            print(f"'{disease}': Neurologic sign match. Score +2. Current score: {scores[disease]}")
            
        # Score for specific cardiac sign (heart murmur)
        if "heart_murmur" in patient_symptoms and "heart_murmur" in profile["symptoms"]:
            scores[disease] += 3
            print(f"'{disease}': Cardiac sign match. Score +3. Current score: {scores[disease]}")

        # Score for geography
        if patient_history["location"] in profile["geography"]:
            scores[disease] += 1
            print(f"'{disease}': Geography match. Score +1. Current score: {scores[disease]}")

        # Score for lab results
        if disease == "C. Borrelia burgdorferi" and patient_labs["Lyme_IgM"] == "positive":
            scores[disease] += 5
            print(f"'{disease}': Positive specific lab. Score +5. Current score: {scores[disease]}")
        
    # Step 4: Identify the most likely diagnosis
    most_likely_diagnosis = max(scores, key=scores.get)

    print("="*30)
    print("Final Scores:")
    for disease, score in scores.items():
        print(f"'{disease}': Final Score = {score}")
    print("="*30)
    
    # Step 5: Print conclusion
    print("Conclusion:")
    print("The patient's presentation with fever, myalgia, neurologic symptoms (disorientation), and cardiac symptoms (heart murmur) after a camping trip is highly characteristic of early disseminated Lyme disease.")
    print("The lab result of an elevated IgM for Lyme confirms an acute infection with Borrelia burgdorferi.")
    print(f"\nThe most likely diagnosis is: {most_likely_diagnosis}")

solve_clinical_case()