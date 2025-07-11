def diagnose_patient():
    """
    Analyzes patient data to determine the most likely diagnosis from a list of choices.
    """
    # Patient's clinical and lab data
    patient_data = {
        "age": 1,
        "symptoms": ["hypertrophic scarring", "erythema", "spasticity"],
        "labs": {"anti-Mi-2": "negative"}
    }

    # Diagnostic knowledge base
    diagnoses = {
        "A. Ectropion": {"score": 0, "reasoning": "Localized eye condition; does not match systemic symptoms."},
        "B. McArdle disease": {"score": 0, "reasoning": "Incorrect age of onset (typically adolescence) and symptoms (lacks skin findings)."},
        "C. Dermatomyositis": {"score": 0, "reasoning": "Matches age, skin, and muscle involvement. Negative anti-Mi-2 is common in juvenile form."},
        "D. McCune Albright syndrome": {"score": 0, "reasoning": "Symptoms do not match the classic triad of this syndrome."},
        "E. Cataracts": {"score": 0, "reasoning": "Localized eye condition; does not match systemic symptoms."}
    }

    # --- Scoring Logic ---
    # This simulates the diagnostic reasoning process.

    # 1. Score for Erythema (skin finding)
    # Dermatomyositis is the only choice where erythema is a key feature.
    score_erythema = 1
    diagnoses["C. Dermatomyositis"]["score"] += score_erythema

    # 2. Score for Muscle Involvement (spasticity)
    # Dermatomyositis is a myopathy (muscle disease).
    score_muscle = 1
    diagnoses["C. Dermatomyositis"]["score"] += score_muscle

    # 3. Score for Age Consistency
    # Juvenile Dermatomyositis fits the patient's age of 1.
    score_age = 1
    diagnoses["C. Dermatomyositis"]["score"] += score_age

    # 4. Score for Lab Result
    # A negative anti-Mi-2 test is common in Juvenile Dermatomyositis.
    score_labs = 1
    diagnoses["C. Dermatomyositis"]["score"] += score_labs

    # --- Print Evaluation ---
    print("Diagnostic Evaluation:")
    print("-" * 30)
    print(f"Patient Profile: Age {patient_data['age']}, Symptoms: {', '.join(patient_data['symptoms'])}, Labs: Anti-Mi-2 Negative")
    print("-" * 30)
    print("Equation for the Most Likely Diagnosis (Dermatomyositis):")
    print(f"Score for Erythema       = {score_erythema}")
    print(f"Score for Muscle Symptom = {score_muscle}")
    print(f"Score for Patient Age    = {score_age}")
    print(f"Score for Lab Finding    = {score_labs}")
    print("------------------------------------------")
    final_score_c = diagnoses['C. Dermatomyositis']['score']
    print(f"Total Score for C          = {final_score_c}")
    
    print("\nConclusion:")
    for diagnosis, data in diagnoses.items():
        if data['score'] > 0:
            print(f"- {diagnosis}: Score {data['score']}. {data['reasoning']}")
        else:
             print(f"- {diagnosis}: Score {data['score']}. {data['reasoning']}")

    most_likely = max(diagnoses, key=lambda k: diagnoses[k]['score'])
    print(f"\nThe most likely diagnosis is: {most_likely}")

# Run the diagnostic analysis
diagnose_patient()