def solve_clinical_case():
    """
    This script analyzes a clinical vignette to determine the most likely
    diagnosis and the corresponding expected imaging finding.
    """

    # Step 1: Define the patient's key symptoms from the vignette.
    patient_symptoms = [
        "transient monocular vision loss (ocular)",
        "pulsatile headaches (neurologic)",
        "joint pain (musculoskeletal)",
        "dyspnea (pulmonary)",
        "hearing loss (neurologic/auditory)",
        "painful lower extremity area (dermatologic)"
    ]

    print("--- Patient's Clinical Presentation ---")
    print("The patient exhibits symptoms across multiple systems:")
    for symptom in patient_symptoms:
        print(f"- {symptom}")

    # Step 2: Determine the most likely diagnosis.
    # The combination of eye, nerve, ear, lung, joint, and skin involvement
    # is highly characteristic of a multisystem granulomatous disease.
    diagnosis = "Sarcoidosis"
    
    print("\n--- Diagnostic Conclusion ---")
    print(f"The constellation of symptoms strongly suggests the diagnosis of: {diagnosis}.")
    print("Specifically, the neurologic symptoms (headaches, vision/hearing loss) point to Neurosarcoidosis.")

    # Step 3: Evaluate the given answer choices.
    answer_choices = {
        'A': 'Periarticular bone demineralization visualized by MRI',
        'B': 'Leptomeningeal enhancement with "snowball" hyperintensities visualized by MRI',
        'C': 'Pleural effusion visualized by chest x-ray',
        'D': 'Vascular hemorrhage visualized by MRI',
        'E': 'Intrasellar mass visualized by MRI'
    }

    print("\n--- Evaluating Answer Choices ---")
    # The key finding for Neurosarcoidosis on an MRI of the brain.
    characteristic_finding = 'Leptomeningeal enhancement with "snowball" hyperintensities visualized by MRI'

    correct_choice_key = None
    for key, description in answer_choices.items():
        if description == characteristic_finding:
            correct_choice_key = key
            break

    if correct_choice_key:
        print(f"The classic imaging finding for Neurosarcoidosis is: '{characteristic_finding}'.")
        print(f"This finding matches answer choice {correct_choice_key}.")
    else:
        print("The characteristic finding was not found in the answer choices.")

if __name__ == "__main__":
    solve_clinical_case()
<<<B>>>