def find_diagnosis():
    """
    Analyzes a clinical vignette to determine the most likely diagnosis.
    """
    # Key findings from the clinical case
    patient_age = 2  # years
    growth_percentile = 10 # %
    key_findings = [
        "pelvic mass",
        "aniridia",
        "elevated blood pressure (hypertension)",
        "delayed speech (developmental delay)"
    ]

    print("Analyzing patient's clinical presentation:")
    print(f"- Patient Age: {patient_age} years old")
    print(f"- Growth: {growth_percentile}th percentile")
    print("- Key Findings:")
    for finding in key_findings:
        print(f"  * {finding}")

    print("\nDiagnostic Reasoning:")
    print("1. The presence of 'aniridia' is a very specific and strong clue.")
    print("2. The combination of aniridia, developmental delay (delayed speech), and an abdominal/pelvic mass is characteristic of WAGR syndrome.")
    print("3. WAGR syndrome is an acronym for Wilms tumor, Aniridia, Genitourinary anomalies, and Retardation.")
    print("4. The 'W' in WAGR, Wilms tumor, is another name for Nephroblastoma.")
    print("5. The patient's hypertension can also be explained by a large kidney tumor (Nephroblastoma) causing increased renin secretion.")
    
    diagnosis = "D. Nephroblastoma"
    
    print("\nConclusion:")
    print(f"Based on the classic association of aniridia with a renal mass in a child, the most likely diagnosis is {diagnosis}.")

find_diagnosis()