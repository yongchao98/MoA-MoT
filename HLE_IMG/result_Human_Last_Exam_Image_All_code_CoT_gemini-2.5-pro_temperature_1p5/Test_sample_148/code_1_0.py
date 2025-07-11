def provide_diagnosis_summary():
    """
    This function summarizes the clinical and pathological findings
    to support the final diagnosis.
    """

    # --- Clinical Information ---
    patient_age = 28
    patient_gender = "Female"
    location = "sella turcica"
    specimen_type = "intraoperative crush smear"

    # --- Key Pathological Findings ---
    key_findings = [
        "Monomorphic cell population arranged in syncytial clusters.",
        "Round to oval nuclei with 'salt-and-pepper' chromatin.",
        "Moderate amount of granular, eosinophilic cytoplasm.",
        "Absence of significant pleomorphism, high mitotic activity, or necrosis."
    ]

    # --- Final Diagnosis ---
    diagnosis = "Pituitary Adenoma"

    # --- Print Report ---
    print("--- Diagnostic Report ---")
    print(f"Patient Information: A {patient_age}-year-old {patient_gender} with a mass in the {location}.")
    print(f"Specimen Type: {specimen_type}")
    print("\nMicroscopic Findings:")
    for finding in key_findings:
        print(f"- {finding}")
    print("\nConclusion:")
    print("The combination of clinical and cytological features is characteristic of a Pituitary Adenoma.")

# Execute the function to display the summary
provide_diagnosis_summary()