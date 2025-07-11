def solve_clinical_case():
    """
    Analyzes a clinical vignette to determine the most likely diagnosis.
    """

    # Step 1: Define patient's clinical presentation
    patient_profile = {
        "History": "53-year-old female with alcohol use disorder and a 40-pack-year smoking history.",
        "Symptoms": "Severe substernal chest pain and pain with swallowing (odynophagia).",
        "Labs": "Elevated C-reactive protein and leukocytosis, indicating inflammation.",
        "Imaging": "Esophageal lumen narrowing and wall thickening.",
        "Endoscopy": "Visually normal mucosa, with no signs of erythema, ulcers, plaques, or strictures."
    }

    print("--- Clinical Case Analysis ---")
    for key, value in patient_profile.items():
        print(f"{key}: {value}")

    print("\n--- Evaluating Differential Diagnoses ---")

    # Step 2: Analyze each diagnosis based on the evidence
    print("\nA. Streptococcal esophagitis & E. Herpes esophagitis:")
    print("   - Reasoning: These are infectious causes. They typically present with visible abnormalities on endoscopy, such as plaques, exudates, or ulcers. The patient's normal endoscopy makes these diagnoses highly unlikely.")

    print("\nD. GERD (Gastroesophageal Reflux Disease):")
    print("   - Reasoning: While GERD can cause chest pain, severe odynophagia and significant wall thickening on imaging are not typical. Furthermore, severe GERD would likely show signs of esophagitis (erythema, erosions) on endoscopy, which are absent here.")
    
    print("\nB. Esophageal adenocarcinoma vs. C. Esophageal squamous cell carcinoma (SCC):")
    print("   - Reasoning: These are the two main types of esophageal cancer. The choice between them often comes down to risk factors.")
    print("   - Adenocarcinoma is most strongly linked to chronic GERD and Barrett's esophagus.")
    print("   - SCC is most strongly linked to smoking and alcohol consumption.")
    print("   - This patient has two MAJOR risk factors for SCC: a 40-pack-year smoking history and alcohol use disorder.")
    print("   - Her symptoms (pain, odynophagia) and imaging findings (wall thickening, narrowing) are highly suspicious for an infiltrating tumor.")

    print("\n--- Reconciling the Endoscopy Result ---")
    print("The key to this case is the 'normal' endoscopy. While it seems reassuring, some esophageal cancers, particularly SCC, can grow within the wall of the esophagus (submucosal or intramural spread) without creating a visible ulcer or mass on the surface initially. The imaging findings of wall thickening are more sensitive for this type of cancer than a visual endoscopic inspection alone.")

    print("\n--- Final Conclusion ---")
    print("Given the overwhelming risk factors (smoking, alcohol) combined with symptoms and imaging findings classic for an infiltrative tumor, Esophageal Squamous Cell Carcinoma is the most likely diagnosis, despite the visually normal endoscopy.")

solve_clinical_case()