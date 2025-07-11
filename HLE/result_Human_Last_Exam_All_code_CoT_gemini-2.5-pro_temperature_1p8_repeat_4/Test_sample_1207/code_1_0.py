import textwrap

def analyze_medical_case():
    """
    Analyzes a clinical vignette to determine the most likely imaging finding.
    """
    
    # Step 1: Deconstruct the patient's clinical presentation.
    patient_profile = {
        "Age": "44 years",
        "Gender": "Female",
        "Primary Symptoms": [
            "Transient monocular vision loss",
            "Pulsatile headaches",
            "Joint pain",
            "Dyspnea (shortness of breath)",
            "Hearing loss"
        ],
        "Developing Symptom": "Painful area on lower extremity"
    }

    print("--- Clinical Presentation Analysis ---")
    print("The patient presents with a multi-system inflammatory condition affecting:")
    print("*   Eyes (monocular vision loss)")
    print("*   Central Nervous System (headaches)")
    print("*   Ears (hearing loss)")
    print("*   Joints (pain)")
    print("*   Lungs (dyspnea)")
    print("*   Skin (painful lower extremity lesion, likely erythema nodosum)")
    print("\nThis constellation of symptoms is highly suggestive of systemic sarcoidosis.")
    print("-" * 35, "\n")

    # Step 2: Evaluate the answer choices based on the likely diagnosis of Sarcoidosis.
    print("--- Evaluation of Imaging Findings ---")
    
    # Choice A
    analysis_a = (
        "A. Periarticular bone demineralization visualized by MRI: This is a finding "
        "in inflammatory arthritis but is not specific to the entire clinical syndrome described."
    )
    print(textwrap.fill(analysis_a, width=80))
    print("-" * 20)
    
    # Choice B
    analysis_b = (
        'B. Leptomeningeal enhancement with "snowball" hyperintensities visualized by MRI: '
        "This is the classic description for neurosarcoidosis. Leptomeningeal enhancement "
        "is caused by inflammation of the brain's lining, leading to headaches. The "
        '"snowball" lesions are sarcoid granulomas in the brain parenchyma. This finding '
        "perfectly aligns with a diagnosis of sarcoidosis, which explains all of the patient's symptoms."
    )
    print(textwrap.fill(analysis_b, width=80))
    print("-" * 20)

    # Choice C
    analysis_c = (
        "C. Pleural effusion visualized by chest x-ray: While possible in pulmonary sarcoidosis, "
        "it's a less common and non-specific finding. Bilateral hilar lymphadenopathy would be more classic."
    )
    print(textwrap.fill(analysis_c, width=80))
    print("-" * 20)

    # Choice D
    analysis_d = (
        "D. Vascular hemorrhage visualized by MRI: This is a possible complication of vasculitis "
        "but is not the primary diagnostic imaging pattern for the underlying disease process, sarcoidosis."
    )
    print(textwrap.fill(analysis_d, width=80))
    print("-" * 20)
    
    # Choice E
    analysis_e = (
        "E. Intrasellar mass visualized by MRI: Sarcoidosis can form a mass in the pituitary region, "
        "but it would not explain the full constellation of symptoms as well as diffuse brain and "
        "meningeal involvement."
    )
    print(textwrap.fill(analysis_e, width=80))
    print("-" * 35, "\n")
    
    # Step 3: Conclude with the best answer.
    print("--- Conclusion ---")
    print("The patient's multi-system symptoms strongly suggest systemic sarcoidosis. The imaging findings "
          "in choice B are pathognomonic for neurosarcoidosis, making it the most expected finding.")

if __name__ == "__main__":
    analyze_medical_case()
    final_answer = "B"
    print(f"\nFinal Answer: <<<B>>>")
