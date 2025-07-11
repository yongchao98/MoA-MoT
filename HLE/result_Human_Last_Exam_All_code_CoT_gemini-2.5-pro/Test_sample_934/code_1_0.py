def solve_medical_case():
    """
    This function encapsulates the reasoning and prints the final answer for the clinical case.
    """

    # Patient Profile
    age = 53
    symptoms = ["substernal chest pain (10/10)", "pain with swallowing"]
    risk_factors = ["alcohol use disorder", "smokes 2 packs/day for over 20 years"]
    labs = ["elevated C-reactive protein", "leukocytosis"]
    imaging = ["esophageal lumen narrowing", "wall thickening"]
    endoscopy = "no signs of erythema, ulcers, plaques, or strictures"

    # Analysis
    # The patient's strong risk factors (heavy smoking, alcohol use) are classic for Esophageal Squamous Cell Carcinoma (SCC).
    # Imaging shows wall thickening, which is highly suggestive of an infiltrative malignancy.
    # While the endoscopy is "normal" on the surface, SCC can grow submucosally (within the wall),
    # which explains the wall thickening on imaging without visible ulcers or masses on the mucosal surface.
    # Other options are less likely:
    # - Infectious esophagitis (A, E) and GERD (D) would typically show mucosal changes (ulcers, erythema) on endoscopy.
    # - Esophageal adenocarcinoma (B) would likely present as a visible mass on endoscopy and is more associated with GERD.
    
    final_diagnosis_choice = "C"
    final_diagnosis_text = "Esophageal squamous cell carcinoma"

    print(f"Based on the analysis, the most likely diagnosis is choice {final_diagnosis_choice}: {final_diagnosis_text}.")

solve_medical_case()