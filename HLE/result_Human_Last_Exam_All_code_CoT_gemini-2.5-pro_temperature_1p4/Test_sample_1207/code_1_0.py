import collections

def find_expected_modality_and_finding():
    """
    This script analyzes a clinical case to determine the most likely imaging finding.
    It breaks down the patient's symptoms, forms a diagnosis, and matches it
    to the most appropriate answer choice.
    """

    # Step 1: Define patient's key symptoms from the case.
    symptoms = {
        "Ocular": "Transient monocular vision loss",
        "Neurologic": ["Pulsatile headaches", "Hearing loss"],
        "Musculoskeletal": "Joint pain",
        "Pulmonary": "Dyspnea",
        "Dermatologic": "Painful lower extremity lesion"
    }

    # Step 2: Analyze the pattern of symptoms.
    print("Clinical Analysis:")
    print("The patient presents with a multisystem inflammatory condition.")
    print(f"Key systems involved are: {', '.join(symptoms.keys())}.")
    print("\nThis combination of neurologic (headache, hearing loss), ocular (vision loss), pulmonary (dyspnea), and musculoskeletal (joint pain) symptoms is highly suggestive of a systemic granulomatous disease, most notably Sarcoidosis.\n")

    # Step 3: Pinpoint the specific diagnosis.
    diagnosis = "Neurosarcoidosis"
    print(f"Given the prominent neurologic and ocular signs, the most likely diagnosis is {diagnosis}.\n")

    # Step 4: Map the diagnosis to expected imaging findings from the choices.
    answer_choices = {
        'A': 'Periarticular bone demineralization visualized by MRI',
        'B': 'Leptomeningeal enhancement with "snowball" hyperintensities visualized by MRI',
        'C': 'Pleural effusion visualized by chest x-ray',
        'D': 'Vascular hemorrhage visualized by MRI',
        'E': 'Intrasellar mass visualized by MRI'
    }

    print("Evaluating Imaging Findings:")
    print(f"The classic and expected MRI finding for {diagnosis} involves inflammation of the brain's coverings (meninges) and the formation of granulomas within the brain tissue.\n")

    # Step 5: Identify the correct answer.
    correct_choice_key = 'B'
    correct_choice_description = answer_choices[correct_choice_key]

    print("Conclusion:")
    print(f"The finding that best matches this pathophysiology is '{correct_choice_description}'.")
    print(f"\nTherefore, the correct answer is Choice {correct_choice_key}.")
    print("--------------------------------------------------")
    print(f"Final Answer: {correct_choice_key}. {correct_choice_description}")


# Execute the diagnostic process
find_expected_modality_and_finding()