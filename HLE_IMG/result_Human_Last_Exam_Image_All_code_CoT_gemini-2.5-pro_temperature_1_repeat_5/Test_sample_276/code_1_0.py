def diagnose_patient():
    """
    Analyzes the patient's clinical and radiological data to determine the most likely diagnosis.
    """
    # Patient Data
    age = 67  # years
    wbc_count = 13000  # cells/mm^3

    # Key Findings
    symptoms = [
        "Chronic (1-month) low-grade fever, weight loss, fatigue, diarrhea",
        "Acute right hemi-abdomen pain and tenderness"
    ]
    history = [
        "Uveitis",
        "Arthritis"
    ]
    lab_findings = [
        f"WBC count of {wbc_count} (leukocytosis)",
        "Positive fecal occult blood test"
    ]
    ct_findings = [
        "Marked circumferential wall thickening of the ileocecal region"
    ]

    # Analysis
    analysis_text = """
The patient's presentation with chronic constitutional symptoms, diarrhea, and acute right-sided abdominal pain,
combined with lab evidence of inflammation (leukocytosis) and bleeding, points to a significant inflammatory
process in the right colon. The CT scan confirms this, showing severe thickening of the ileocecal region.

While several conditions can cause this (e.g., infectious colitis, TB, lymphoma), the patient's
past medical history of both uveitis and arthritis is highly specific for an extra-intestinal
manifestation of Inflammatory Bowel Disease (IBD). Given that Crohn's Disease classically affects
the ileocecal region, it provides the best unifying diagnosis for all the patient's findings.
"""

    # Final Diagnosis
    diagnosis = "A. Crohn's Disease"

    print("--- Patient Case Analysis ---")
    print("\nKey Symptoms:")
    for s in symptoms:
        print(f"- {s}")
    
    print("\nRelevant History:")
    for h in history:
        print(f"- {h}")

    print("\nKey Lab and Exam Findings:")
    for f in lab_findings:
        print(f"- {f}")

    print("\nKey CT Findings:")
    for f in ct_findings:
        print(f"- {f}")
    
    print("\n--- Diagnostic Reasoning ---")
    print(analysis_text)
    
    print("--- Conclusion ---")
    print(f"The most likely diagnosis is: {diagnosis}")

if __name__ == "__main__":
    diagnose_patient()