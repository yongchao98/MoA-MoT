import textwrap

def solve_clinical_case():
    """
    This function analyzes a clinical vignette and determines the most likely imaging finding.
    """

    # Step 1: Analyze the key symptoms from the case description.
    patient_profile = "44-year-old woman"
    symptoms = [
        "Transient monocular vision loss",
        "Pulsatile headaches",
        "Joint pain",
        "Dyspnea (shortness of breath)",
        "Hearing loss",
        "Painful lower extremity area"
    ]

    # Step 2: Explain the diagnostic reasoning.
    reasoning_text = """
The patient presents with a classic multisystem inflammatory disorder affecting the eyes, nervous system, lungs, joints, and skin. This constellation of symptoms strongly points towards a systemic granulomatous disease.

The most likely unifying diagnosis is Sarcoidosis. Here's how it explains each symptom:
- Neurosarcoidosis (Sarcoidosis of the nervous system) can cause headaches, inflammation of the optic nerve (leading to vision loss), and involvement of the auditory nerve (causing hearing loss).
- Pulmonary Sarcoidosis, the most common form, explains the dyspnea.
- Sarcoid Arthropathy explains the joint pain.
- Erythema Nodosum, a common skin manifestation of sarcoidosis, presents as painful nodules on the lower extremities.

Given the strong evidence for Neurosarcoidosis, we must identify its characteristic imaging finding.
"""

    # Step 3: Evaluate the provided answer choices.
    evaluation_text = """
Evaluating the Imaging Modalities:
- A. Periarticular bone demineralization is characteristic of rheumatoid arthritis, not sarcoidosis.
- B. Leptomeningeal enhancement (inflammation of the brain's lining) and parenchymal 'snowball' hyperintensities on MRI are classic and highly specific findings for Neurosarcoidosis. This matches our diagnosis perfectly.
- C. Pleural effusion is a non-specific finding and less common in sarcoidosis than interstitial lung disease or hilar adenopathy. It doesn't explain the neurological symptoms.
- D. Vascular hemorrhage is also non-specific and does not point directly to the underlying systemic disease as well as choice B does.
- E. An intrasellar (pituitary) mass would not account for the patient's widespread systemic symptoms (dyspnea, joint pain, hearing loss, skin lesion).
"""

    # Step 4: Conclude and state the final answer.
    conclusion_text = "The most fitting diagnosis is Sarcoidosis, and the expected imaging finding for its neurological involvement is described in choice B."

    # Print the full analysis.
    print("Clinical Reasoning:\n" + "="*20)
    print(textwrap.dedent(reasoning_text))
    print("Evaluation of Answer Choices:\n" + "="*30)
    print(textwrap.dedent(evaluation_text))
    print("Conclusion:\n" + "="*12)
    print(conclusion_text)
    print("\nTherefore, the correct answer is:")
    print("B. Leptomeningeal enhancement with 'snowball' hyperintensities visualized by MRI")

# Execute the function to print the solution.
solve_clinical_case()