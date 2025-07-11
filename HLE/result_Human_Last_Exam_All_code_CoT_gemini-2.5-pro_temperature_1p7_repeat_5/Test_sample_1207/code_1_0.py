import textwrap

def solve_clinical_case():
    """
    Analyzes a clinical vignette to determine the most likely imaging finding.
    """
    
    print("Analyzing the patient's case step-by-step:\n")

    # Step 1: Deconstruct the patient's symptoms
    step1_title = "Step 1: Symptom Analysis"
    step1_text = """
    The patient presents with a constellation of symptoms affecting multiple organ systems:
    - Ocular: Transient monocular vision loss.
    - Neurological: Pulsatile headaches, hearing loss.
    - Rheumatological: Joint pain.
    - Pulmonary: Dyspnea (shortness of breath).
    - Dermatological: Painful area on the lower extremity (suggestive of a lesion like erythema nodosum).
    This multi-system involvement points towards a systemic inflammatory or autoimmune disease.
    """
    print(f"--- {step1_title} ---\n{textwrap.dedent(step1_text)}")

    # Step 2: Form a differential diagnosis
    step2_title = "Step 2: Differential Diagnosis"
    step2_text = """
    Given the symptoms, a strong candidate is Sarcoidosis, a multi-system inflammatory disease characterized by the formation of granulomas.
    - Neurosarcoidosis can affect the cranial nerves (optic and vestibulocochlear nerves), causing the vision and hearing loss, as well as inflammation of the meninges (leptomeninges), causing headaches.
    - Pulmonary sarcoidosis is a classic cause of dyspnea.
    - Sarcoidosis is also known to cause arthritis (joint pain) and skin lesions.
    """
    print(f"--- {step2_title} ---\n{textwrap.dedent(step2_text)}")

    # Step 3: Evaluate the imaging options
    step3_title = "Step 3: Evaluating the Imaging Findings"
    step3_text = """
    We will now evaluate each choice based on its fit with a diagnosis of systemic sarcoidosis:
    - A. Periarticular bone demineralization: Non-specific finding of joint inflammation, doesn't explain the main neurological issues.
    - C. Pleural effusion: Can be caused by sarcoidosis, but is less specific and doesn't explain the neurological symptoms.
    - D. Vascular hemorrhage: This is a potential complication, not the primary imaging sign of the underlying chronic inflammatory disease.
    - E. Intrasellar mass: Would not cause monocular vision loss or the other systemic symptoms.
    
    - B. Leptomeningeal enhancement with 'snowball' hyperintensities visualized by MRI: This finding is highly characteristic of neurosarcoidosis. The leptomeningeal enhancement reflects inflammation around the brain and spinal cord, causing headaches. Granulomatous inflammation can affect the cranial nerves. 'Snowball' hyperintensities represent parenchymal granulomas in the brain. This finding provides a single, unifying diagnosis for the patient's most specific symptoms.
    """
    print(f"--- {step3_title} ---\n{textwrap.dedent(step3_text)}")
    
    # Final conclusion
    final_conclusion_title = "Final Conclusion"
    final_conclusion_text = """
    The patient's multi-system symptoms are best explained by Sarcoidosis. The most specific and expected imaging finding for Neurosarcoidosis among the given choices is B.
    """
    print(f"--- {final_conclusion_title} ---\n{textwrap.dedent(final_conclusion_text)}")

solve_clinical_case()
print("<<<B>>>")