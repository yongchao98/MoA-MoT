def diagnose_patient():
    """
    This function analyzes the patient's clinical and radiological data to arrive at the most likely diagnosis.
    """
    # Patient Profile
    age = 67
    sex = "Female"
    
    # Key Symptoms (chronic, >1 month)
    symptoms = ["low-grade fever", "weight loss", "fatigue", "diarrhea"]
    
    # Acute Symptom
    acute_symptom = "acute right-sided abdominal pain"
    
    # Significant Past Medical History
    past_history = ["uveitis", "arthritis"]
    
    # Exam and Lab Findings
    wbc_count = 13000  # cells/mcL
    fecal_occult_blood = "Positive"
    exam_findings = "Right hemiabdomen tenderness and guarding"
    
    # CT Scan Findings
    ct_findings = "Marked wall thickening of the ileocecal region"
    
    # Analysis
    print("Analyzing the case:")
    print(f"A {age}-year-old {sex} presents with a month-long history of constitutional symptoms and diarrhea.")
    print(f"Key extraintestinal manifestations in her history include: {', '.join(past_history)}.")
    print(f"Labs show inflammation (WBC: {wbc_count}) and GI bleeding (Fecal Occult Blood: {fecal_occult_blood}).")
    print(f"The physical exam and CT scan localize the pathology to the ileocecal region, showing {ct_findings}.")
    
    # Differential Diagnosis Logic
    print("\nEvaluating top differentials:")
    print(" - Ileocecal Tuberculosis: Can mimic symptoms and CT findings, but doesn't typically explain the history of uveitis and arthritis.")
    print(" - GI Lymphoma: Can cause wall thickening and constitutional symptoms, but the autoimmune-related history (uveitis, arthritis) is less typical.")
    print(" - Crohn's Disease: Classically involves the ileocecal region, causes chronic inflammation, weight loss, diarrhea, and abdominal pain. Crucially, uveitis and arthritis are well-known extraintestinal manifestations.")
    
    # Conclusion
    final_diagnosis = "Crohn's Disease"
    print(f"\nConclusion: The combination of chronic ileocecal inflammation with a history of uveitis and arthritis makes {final_diagnosis} the most likely diagnosis.")

diagnose_patient()