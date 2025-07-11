def diagnose_patient():
    """
    Analyzes patient symptoms to determine the most likely diagnosis.
    """
    patient = {
        "age": "1-year-old",
        "symptoms": ["hypertrophic scarring", "erythema", "spasticity"],
        "labs": {"anti-Mi-2": "negative"}
    }

    diagnoses = {
        "A. Ectropion": "Isolated eye condition (eyelid turning outward). Does not explain skin or muscle symptoms.",
        "B. McArdle disease": "Metabolic muscle disease, presents with exercise intolerance. Does not explain the skin findings.",
        "C. Dermatomyositis": "Inflammatory disease of skin and muscle. Explains erythema (skin) and spasticity (muscle involvement). Severe skin manifestations can lead to scarring. A negative anti-Mi-2 test is common in the juvenile form.",
        "D. McCune Albright syndrome": "Genetic disorder affecting bones, skin (café-au-lait spots), and endocrine system. Does not match the presented symptoms.",
        "E. Cataracts": "Isolated eye condition (clouding of the lens). Does not explain skin or muscle symptoms."
    }

    print("Analyzing patient case:")
    print(f" - Patient Age: {patient['age']}")
    print(f" - Symptoms: {', '.join(patient['symptoms'])}")
    print(f" - Lab Results: Anti-Mi-2 is {patient['labs']['anti-Mi-2']}\n")

    print("Evaluating potential diagnoses:")
    best_match = None
    highest_score = 0

    for diagnosis, explanation in diagnoses.items():
        score = 0
        if "skin" in explanation and "muscle" in explanation:
            # This disease affects both systems presented by the patient
            score += 2
        if "erythema" in explanation.lower() and "erythema" in patient["symptoms"]:
            score += 1
        if "negative anti-Mi-2" in explanation.lower() and patient["labs"]["anti-Mi-2"] == "negative":
            # This specific lab finding is consistent with this diagnosis
            score += 1

        print(f" - {diagnosis}: {explanation}")
        if score > highest_score:
            highest_score = score
            best_match = diagnosis

    print("\nConclusion:")
    print(f"The most likely diagnosis is {best_match} because it provides a unified explanation for both the skin (erythema, scarring) and muscle (spasticity) findings in a pediatric patient. The negative anti-Mi-2 lab is also consistent with this diagnosis.")

diagnose_patient()