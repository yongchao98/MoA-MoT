def diagnose_esophageal_condition():
    """
    Analyzes patient data to determine the most likely diagnosis for esophageal symptoms.
    """

    # --- Patient Profile & Clinical Findings ---
    patient_age = 53
    smoking_packs_per_day = 2
    smoking_years = 20
    
    # Key Factors for Diagnosis
    risk_factors = ["heavy_smoking", "alcohol_use_disorder"]
    symptoms = ["substernal_chest_pain", "odynophagia"] # odynophagia is pain with swallowing
    imaging_findings = ["esophageal_wall_thickening", "lumen_narrowing"]
    endoscopy_findings = "normal_mucosa" # No ulcers, plaques, or erythema

    print("Analyzing the clinical case of a {}-year-old woman.".format(patient_age))
    print("Significant Risk Factors: {} and {}.".format(risk_factors[0], risk_factors[1]))
    print("Presenting Symptoms: {} and {}.".format(symptoms[0], symptoms[1]))
    print("Key Test Results:")
    print("  - Imaging: {} and {}.".format(imaging_findings[0], imaging_findings[1]))
    print("  - Endoscopy: {}.".format(endoscopy_findings))
    print("-" * 20)
    print("Evaluating the differential diagnosis:")

    # --- Diagnostic Logic ---
    
    # Evaluate infectious and inflammatory causes
    if endoscopy_findings == "normal_mucosa":
        print("-> [Rule Out A, D, E]: The normal endoscopic findings make infectious causes like Streptococcal (A) or Herpes (E) esophagitis highly unlikely, as they typically present with ulcers or plaques.")
        print("-> Similarly, severe GERD (D) causing these symptoms would likely show esophagitis (erythema, erosions) on endoscopy.")
    
    # Differentiate between cancers
    print("-> Both Esophageal Squamous Cell Carcinoma (SCC) and Adenocarcinoma are considerations.")
    print("-> Patient has major risk factors for SCC: smoking {} packs/day for {} years and alcohol use.".format(smoking_packs_per_day, smoking_years))
    
    # The deciding factor
    if "esophageal_wall_thickening" in imaging_findings and endoscopy_findings == "normal_mucosa":
        print("\n-> The critical finding is wall thickening on imaging despite a normal-appearing mucosa on endoscopy.")
        print("-> This pattern is characteristic of an infiltrative tumor, one that grows within the esophageal wall without breaking through the surface mucosa.")
        print("-> This presentation is strongly associated with Esophageal Squamous Cell Carcinoma (C).")
        print("-" * 20)
        final_diagnosis = "C. Esophageal squamous cell carcinoma"
    else:
        final_diagnosis = "Inconclusive based on the provided logic."
        
    print("Most Likely Diagnosis:", final_diagnosis)

diagnose_esophageal_condition()
<<<C>>>