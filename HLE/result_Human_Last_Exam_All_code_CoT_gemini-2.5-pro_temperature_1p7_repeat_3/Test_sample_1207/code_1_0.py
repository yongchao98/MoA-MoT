def solve_clinical_case():
    """
    Analyzes a clinical vignette to determine the most likely diagnosis and corresponding imaging finding.
    """
    # 1. Define patient symptoms from the case description.
    symptoms = [
        "transient monocular vision loss",
        "pulsatile headaches",
        "joint pain",
        "dyspnea",
        "hearing loss",
        "painful lower extremity lesion"
    ]

    # 2. Define potential diagnoses and their corresponding imaging findings from the answer choices.
    diagnoses_and_findings = {
        "Rheumatoid Arthritis": "A. Periarticular bone demineralization visualized by MRI",
        "Neurosarcoidosis": "B. Leptomeningeal enhancement with \"snowball\" hyperintensities visualized by MRI",
        "Pulmonary Disease (non-specific)": "C. Pleural effusion visualized by chest x-ray",
        "Cerebral Hemorrhage": "D. Vascular hemorrhage visualized by MRI",
        "Pituitary Adenoma": "E. Intrasellar mass visualized by MRI"
    }

    # 3. Analyze the symptoms to find the best-fit diagnosis.
    # The combination of neurological, ocular, auditory, pulmonary, cutaneous, and joint symptoms
    # in a middle-aged woman is highly characteristic of a multisystem granulomatous disease.
    best_fit_diagnosis = "Sarcoidosis (presenting as Neurosarcoidosis)"

    # 4. Print the step-by-step reasoning.
    print("--- Clinical Case Analysis ---")
    print(f"Patient Symptoms: {', '.join(symptoms)}.")
    print("\nReasoning:")
    print("The patient presents with a multisystem inflammatory condition. The constellation of symptoms strongly points towards Sarcoidosis.")
    print("- Ocular (vision loss), Neurological (headaches), and Auditory (hearing loss) symptoms suggest Neurosarcoidosis.")
    print("- Dyspnea points to pulmonary involvement (classic for Sarcoidosis).")
    print("- Joint pain and painful skin lesions are also common systemic features of Sarcoidosis.")

    print("\nEvaluating the Imaging Findings:")
    # Find the imaging finding that matches the best-fit diagnosis.
    for diagnosis, finding in diagnoses_and_findings.items():
        if "Neurosarcoidosis" in diagnosis:
            correct_finding = finding
            print(f"- {finding}: This is a classic finding for {diagnosis}, which aligns perfectly with the clinical picture.")
        else:
            print(f"- {finding}: This finding corresponds to {diagnosis}, which is less likely to explain the full range of the patient's symptoms.")

    # 5. Output the final answer.
    print("\n--- Conclusion ---")
    print("The most likely diagnosis is Sarcoidosis, and the expected imaging finding is the one characteristic of its neurological manifestation.")
    final_answer_choice = correct_finding.split('.')[0]
    print(f"\nThe correct imaging modality and finding is: {correct_finding}")
    print(f"Final Answer Choice: {final_answer_choice}")


solve_clinical_case()