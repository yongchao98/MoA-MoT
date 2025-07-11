def solve_diagnosis():
    """
    This function analyzes the clinical case to determine the most likely diagnosis.
    """
    # Patient Information
    age = 67
    sex = "female"

    # Symptoms
    symptoms = {
        "duration_months": 1,
        "fever": "low-grade",
        "weight_loss": True,
        "fatigue": True,
        "diarrhea": True,
        "acute_pain": "right hemiabdomen"
    }

    # Past Medical History (PMH) - Key Clue
    pmh = ["gallstones", "uveitis", "arthritis"]

    # Lab and Exam Findings
    wbc_count = 13000  # cells/mm^3
    fecal_occult_blood = "positive"
    exam_findings = "tenderness and guarding at right hemiabdomen"

    # CT Findings
    ct_findings = "marked wall thickening of the cecum and terminal ileum"

    # Analysis
    # The key to this diagnosis is linking the gastrointestinal presentation
    # with the patient's past medical history.
    # The combination of ileitis (inflammation of the ileum, suggested by RLQ pain and CT)
    # with extraintestinal manifestations like uveitis and arthritis is highly characteristic
    # of Crohn's Disease.

    # Differential Diagnosis Consideration:
    # - Ileocecal Tuberculosis (C) is a strong mimic, but the presence of both uveitis and
    #   arthritis makes Crohn's more probable.
    # - Other infectious causes (B, D, E) are less likely given the chronicity and PMH.
    # - Structural/Vascular issues (F, G, H, I) do not fit the overall clinical picture.
    # - Neoplasm (K) is possible but less likely than Crohn's given the specific PMH.

    final_diagnosis_code = "A"
    final_diagnosis_name = "Crohn's Disease"

    print(f"Patient Profile: {age}-year-old {sex}")
    print(f"Key Symptoms: Chronic diarrhea, weight loss, fatigue for {symptoms['duration_months']} month, with acute right-sided abdominal pain.")
    print(f"Key Exam/Lab Findings: WBC {wbc_count}, Positive Fecal Occult Blood, Right Hemiabdomen Tenderness.")
    print(f"Key Past Medical History: {', '.join(pmh)}")
    print(f"Key CT Finding: {ct_findings}")
    print("\nConclusion:")
    print("The constellation of findings, particularly the involvement of the ileocecal region")
    print("combined with a classic history of extraintestinal manifestations (uveitis and arthritis),")
    print(f"makes '{final_diagnosis_name}' the most likely diagnosis.")
    print(f"The correct answer choice is {final_diagnosis_code}.")

solve_diagnosis()