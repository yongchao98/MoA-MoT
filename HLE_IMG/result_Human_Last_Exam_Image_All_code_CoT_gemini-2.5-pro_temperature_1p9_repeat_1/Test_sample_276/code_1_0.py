def solve_medical_case():
    """
    This function analyzes the provided medical case and determines the most likely diagnosis.
    """

    # Clinical Information Summary
    patient_age = 67
    symptoms = [
        "low-grade fever",
        "weight loss",
        "fatigue",
        "diarrhea (1 month)",
        "acute right hemiabdomen pain"
    ]
    past_medical_history = ["gallstones", "uveitis", "arthritis"]
    exam_findings = ["right hemiabdomen tenderness", "guarding"]
    lab_results = {
        "WBC Count": 13000, # Normal is 4,500-11,000
        "Fecal Occult Blood Test": "Positive"
    }
    ct_findings = [
        "marked wall thickening of terminal ileum and cecum (ileocecal region)",
        "surrounding mesenteric inflammation"
    ]

    # Diagnosis Rationale:
    # The constellation of findings (chronic GI and constitutional symptoms, right lower quadrant focus,
    # and CT showing ileocecal inflammation) points towards an inflammatory process.
    # The key to distinguishing between mimics like infection (Yersinia, TB) and other conditions
    # is the patient's past medical history.
    # Uveitis and arthritis are well-known extra-intestinal manifestations of Inflammatory Bowel Disease.
    # Crohn's Disease is the form of IBD that classically involves the ileocecal region.
    # Therefore, Crohn's Disease is the most unifying diagnosis.

    diagnosis = "A. Crohn's Disease"
    
    print("Patient Profile:")
    print(f"Age: {patient_age} year old female")
    print(f"Symptoms: {', '.join(symptoms)}")
    print(f"Key Past Medical History: {', '.join(past_medical_history)}")
    print(f"Key Lab/Exam Findings: WBC {lab_results['WBC Count']}, Positive Fecal Occult Blood, {', '.join(exam_findings)}")
    print(f"Key CT Findings: {', '.join(ct_findings)}")
    print("\nAnalysis:")
    print("The patient's clinical presentation with chronic diarrhea, weight loss, fever, and acute-on-chronic right lower quadrant pain, combined with CT findings of inflammation in the ileocecal region, is highly suggestive of inflammatory bowel disease.")
    print("The patient's history of both uveitis and arthritis provides strong evidence for a systemic inflammatory condition. These are classic extra-intestinal manifestations of Crohn's Disease.")
    print("While infections (like Yersinia or TB) or malignancy can mimic this presentation, the presence of these specific extra-intestinal manifestations makes Crohn's Disease the most likely diagnosis.")
    
    # Final answer choice letter
    final_answer_letter = diagnosis.split('.')[0]
    print(f"\nMost Likely Diagnosis: {diagnosis}")
    
solve_medical_case()