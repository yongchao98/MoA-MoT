def diagnose_patient():
    """
    This function simulates the diagnostic process for the given patient case
    by scoring potential diagnoses based on clinical and imaging findings.
    """

    # --- Patient Data ---
    patient_findings = {
        "age": 67,
        "duration": "chronic",  # >4 weeks
        "symptoms": {"low-grade fever", "weight loss", "fatigue", "diarrhea", "abdominal pain"},
        "history": {"gallstones", "uveitis", "arthritis", "animal_contact"},
        "exam": {"right_hemiabdomen_tenderness", "guarding"},
        "labs": {"leukocytosis", "fecal_occult_blood_positive"},
        "imaging": {"ileocecal_thickening"}
    }

    # --- Medical Knowledge Base ---
    # We define features for each diagnosis. A strong positive feature gets +2, a typical one +1.
    # A strong negative feature (contradiction) gets -5.
    diagnoses_kb = {
        "A. Crohn's Disease": {"score": 0, "features": {
            "duration_chronic": 2, "symptoms_constitutional": 2, "symptoms_gi": 1,
            "history_uveitis": 5, "history_arthritis": 5, # Very high weight
            "exam_rlq": 2, "labs_inflammatory": 1, "labs_bleeding": 1, "imaging_ileocecal": 2}},
        "B. Yersinia Colitis": {"score": 0, "features": {
            "duration_chronic": -5, "symptoms_gi": 1, "exam_rlq": 1, "history_animal_contact": 1, "imaging_ileocecal": 2}},
        "C. Ileocecal Tuberculosis": {"score": 0, "features": {
            "duration_chronic": 2, "symptoms_constitutional": 2, "symptoms_gi": 1,
            "exam_rlq": 2, "labs_inflammatory": 1, "labs_bleeding": 1, "imaging_ileocecal": 2}},
        "D. Salmonella Enteritis": {"score": 0, "features": {"duration_chronic": -5, "symptoms_gi": 1}},
        "E. C. difficile Colitis": {"score": 0, "features": {"duration_chronic": -5, "symptoms_gi": 1, "labs_inflammatory": 1}},
        "F. Cecal Volvulus": {"score": 0, "features": {"duration_chronic": -5, "imaging_ileocecal": -5}}, # Image is not volvulus
        "G. Ischemic Colitis": {"score": 0, "features": {"duration_chronic": -5, "symptoms_gi": 1}},
        "H. Pneumoperitoneum": {"score": 0, "features": {"imaging_ileocecal": -5}}, # Not a diagnosis of cause
        "I. Ovarian Torsion": {"score": 0, "features": {"duration_chronic": -5}},
        "J. Celiac Disease": {"score": 0, "features": {"imaging_ileocecal": -5, "symptoms_constitutional": 1}},
        "K. Gastrointestinal Lymphoma": {"score": 0, "features": {
            "duration_chronic": 2, "symptoms_constitutional": 2, "symptoms_gi": 1,
            "exam_rlq": 1, "labs_bleeding": 1, "imaging_ileocecal": 2}}
    }

    print("Analyzing patient findings to determine the most likely diagnosis...\n")

    # --- Scoring Logic ---
    print("Scoring diagnoses based on key findings:")

    # Duration
    if patient_findings["duration"] == "chronic":
        diagnoses_kb["A. Crohn's Disease"]["score"] += diagnoses_kb["A. Crohn's Disease"]["features"]["duration_chronic"]
        diagnoses_kb["B. Yersinia Colitis"]["score"] += diagnoses_kb["B. Yersinia Colitis"]["features"]["duration_chronic"]
        diagnoses_kb["C. Ileocecal Tuberculosis"]["score"] += diagnoses_kb["C. Ileocecal Tuberculosis"]["features"]["duration_chronic"]
        diagnoses_kb["K. Gastrointestinal Lymphoma"]["score"] += diagnoses_kb["K. Gastrointestinal Lymphoma"]["features"]["duration_chronic"]
        print(f"- Chronic presentation (+{diagnoses_kb['A. Crohn\'s Disease']['features']['duration_chronic']} for Crohn's, TB, Lymphoma; {diagnoses_kb['B. Yersinia Colitis']['features']['duration_chronic']} for acute infections)")

    # Location and Imaging
    if "ileocecal_thickening" in patient_findings["imaging"]:
        for dx in ["A. Crohn's Disease", "B. Yersinia Colitis", "C. Ileocecal Tuberculosis", "K. Gastrointestinal Lymphoma"]:
             diagnoses_kb[dx]["score"] += diagnoses_kb[dx]["features"]["imaging_ileocecal"]
        print(f"- Ileocecal thickening on CT (+{diagnoses_kb['A. Crohn\'s Disease']['features']['imaging_ileocecal']} points for Crohn's, Yersinia, TB, Lymphoma)")

    # Constitutional Symptoms
    if {"low-grade fever", "weight loss"}.issubset(patient_findings["symptoms"]):
        for dx in ["A. Crohn's Disease", "C. Ileocecal Tuberculosis", "K. Gastrointestinal Lymphoma"]:
            diagnoses_kb[dx]["score"] += diagnoses_kb[dx]["features"]["symptoms_constitutional"]
        print(f"- Constitutional symptoms (+{diagnoses_kb['A. Crohn\'s Disease']['features']['symptoms_constitutional']} points for Crohn's, TB, Lymphoma)")

    # The CRITICAL finding: Extra-intestinal manifestations
    if "uveitis" in patient_findings["history"] and "arthritis" in patient_findings["history"]:
        crohn_bonus = diagnoses_kb["A. Crohn's Disease"]["features"]["history_uveitis"] + diagnoses_kb["A. Crohn's Disease"]["features"]["history_arthritis"]
        diagnoses_kb["A. Crohn's Disease"]["score"] += crohn_bonus
        print(f"- History of Uveitis AND Arthritis (CRITICAL FINDING): +{crohn_bonus} points for Crohn's Disease.")

    # --- Final Calculation ---
    print("\n--- Final Scores ---")
    most_likely_diagnosis = ""
    max_score = -99

    for diagnosis, data in sorted(diagnoses_kb.items(), key=lambda item: item[1]['score'], reverse=True):
        score = data['score']
        print(f"{diagnosis}: {score}")
        if score > max_score:
            max_score = score
            most_likely_diagnosis = diagnosis

    print("\n--- Conclusion ---")
    print(f"The analysis indicates that the most likely diagnosis is {most_likely_diagnosis.split('.')[1].strip()}.")
    print("This is due to the perfect alignment of chronic ileocecal inflammation with the highly specific extra-intestinal manifestations (uveitis and arthritis), which are classic for this condition.")

    # Return the final answer in the specified format
    final_answer_letter = most_likely_diagnosis.split('.')[0]
    return final_answer_letter

if __name__ == '__main__':
    answer = diagnose_patient()
    print(f"<<<{answer}>>>")