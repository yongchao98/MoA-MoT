def analyze_clinical_case():
    """
    Analyzes a clinical vignette to determine the most likely diagnosis.
    """
    # Patient Data
    patient_age = 1  # years
    key_findings = ["hypertrophic scarring", "erythema", "spasticity"]
    lab_results = {"anti-Mi-2": "negative"}

    # Diagnosis Information
    diagnoses = {
        "A. Ectropion": "Localized eyelid condition. Does not match systemic symptoms.",
        "B. McArdle disease": "Metabolic muscle disease, different presentation and onset age.",
        "C. Dermatomyositis": "Inflammatory skin/muscle disease. Juvenile form fits age, skin signs (erythema, scarring from ulceration), and muscle signs (spasticity). Negative anti-Mi-2 is common in juveniles.",
        "D. McCune Albright syndrome": "Specific triad of bone, skin spots, and endocrine issues not seen here.",
        "E. Cataracts": "Localized eye condition. Does not match systemic symptoms."
    }

    print("Patient Profile:")
    print(f"Age: {patient_age} year old")
    print(f"Key Findings: {', '.join(key_findings)}")
    print(f"Lab Result: anti-Mi-2 is {lab_results['anti-Mi-2']}")
    print("-" * 30)

    print("Analysis:")
    print("The patient's symptoms involve multiple systems (skin and musculoskeletal), suggesting a systemic disease.")
    print("Choices A and E are localized conditions and can be ruled out.")
    print("Choices B and D are systemic but their classic presentations do not match the patient's findings.")
    print("\nChoice C, Dermatomyositis (specifically the juvenile form), accounts for all findings:")
    print("  - Erythema (skin inflammation)")
    print("  - Hypertrophic scarring (can result from skin vasculopathy/ulceration)")
    print("  - Spasticity (muscle involvement)")
    print("  - The patient's age of 1 is consistent with Juvenile Dermatomyositis.")
    print("  - A negative anti-Mi-2 test does not rule out the diagnosis in children.")

    most_likely = "C. Dermatomyositis"
    print("\nConclusion: The most likely diagnosis is " + most_likely)

analyze_clinical_case()