def diagnose_tick_borne_illness():
    """
    Analyzes a clinical vignette to determine the most likely diagnosis.
    """
    patient_data = {
        "location": "Oklahoma",
        "symptoms": ["fever", "headaches", "myalgia", "disorientation", "heart murmur"],
        "labs": "Positive IgM Lyme, Negative IgG Lyme"
    }

    diseases = {
        "A. Babesia microti": {
            "common_geography": "Northeast and Upper Midwest US",
            "match_score": 2,
            "notes": "Geography is a poor match."
        },
        "B. Plasmodium": {
            "common_geography": "Malaria-endemic zones (not Oklahoma)",
            "match_score": 0,
            "notes": "Incorrect vector and geography."
        },
        "C. Borrelia burgdorferi (Lyme)": {
            "common_geography": "Northeast and Upper Midwest US",
            "match_score": 3,
            "notes": "Symptoms fit, but geography is less likely. The positive IgM could be the true result or a cross-reaction."
        },
        "D. Ehrlichia": {
            "common_geography": "South-Central and Eastern US (including Oklahoma)",
            "match_score": 5,
            "notes": "Excellent match for geography and symptoms. Known to cross-react with Lyme serology, explaining the lab results."
        },
        "E. Rickettsia rickettsii (RMSF)": {
            "common_geography": "Throughout US, including Oklahoma",
            "match_score": 4,
            "notes": "Excellent match for geography and symptoms. A very strong possibility, but the Lyme cross-reactivity points more specifically to Ehrlichia."
        }
    }

    print("Analyzing the clinical case:")
    print(f"Patient Location: {patient_data['location']}")
    print(f"Patient Symptoms: {', '.join(patient_data['symptoms'])}")
    print(f"Patient Labs: {patient_data['labs']}\n")

    print("Evaluating potential diagnoses:")
    best_match = ""
    highest_score = -1

    for disease, data in diseases.items():
        print(f"- {disease}:")
        print(f"  Geography: {data['common_geography']}")
        print(f"  Notes: {data['notes']}")
        if data['match_score'] > highest_score:
            highest_score = data['match_score']
            best_match = disease

    print("\nConclusion:")
    print("The patient's presentation in Oklahoma with fever, headache, myalgia, and disorientation is classic for a tick-borne illness common to that region.")
    print("The positive Lyme IgM is likely a serologic cross-reaction, a known phenomenon with Ehrlichiosis.")
    print(f"Therefore, the most likely diagnosis is: {best_match}")

diagnose_tick_borne_illness()