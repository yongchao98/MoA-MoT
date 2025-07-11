def solve_clinical_case():
    """
    Analyzes the clinical vignette to determine the most likely diagnosis and corresponding imaging finding.
    """
    patient_symptoms = {
        "Ocular": "Transient monocular vision loss",
        "Neurologic": "Pulsatile headaches, hearing loss",
        "Musculoskeletal": "Joint pain",
        "Pulmonary": "Dyspnea",
        "Cutaneous": "Painful area on lower extremity"
    }

    # The combination of multisystem symptoms (eyes, nerves, lungs, joints, skin)
    # strongly suggests a systemic granulomatous disease like Sarcoidosis.
    # Neurosarcoidosis specifically explains the headache, vision, and hearing loss.
    # The painful leg lesion is likely erythema nodosum, also associated with sarcoidosis.

    reasoning = """
The patient's presentation with transient monocular vision loss, headaches, hearing loss, joint pain, dyspnea, and a painful lower extremity lesion is highly suggestive of a multisystem inflammatory disorder. Sarcoidosis is a prime candidate.
- Neurosarcoidosis can cause cranial nerve palsies (affecting vision and hearing) and meningeal inflammation (causing headaches).
- Pulmonary sarcoidosis causes dyspnea.
- Joint and skin involvement are also common.
The imaging finding most specific for neurosarcoidosis among the choices is leptomeningeal enhancement with periventricular "snowball" hyperintensities on an MRI of the brain.
"""

    answer_choice = "B"
    answer_description = "Leptomeningeal enhancement with \"snowball\" hyperintensities visualized by MRI"

    print("--- Clinical Reasoning ---")
    print(reasoning)
    print("--- Conclusion ---")
    print(f"The most expected image modality and finding is: {answer_choice}. {answer_description}")
    print("\nFinal Answer Code:")
    print(f'<<<{answer_choice}>>>')

solve_clinical_case()