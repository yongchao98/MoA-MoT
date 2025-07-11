def diagnose_esophageal_condition():
    """
    Analyzes a clinical vignette to determine the most likely diagnosis.
    The code will print the reasoning step-by-step.
    """
    
    # Key data points from the case
    age = 53
    smoking_packs_per_day = 2
    smoking_years = 20
    pack_years = smoking_packs_per_day * smoking_years

    print("Step 1: Evaluating Patient's Risk Factors")
    print(f" - The patient is a {age}-year-old woman.")
    print(" - Major risk factor: Heavy smoking.")
    print(f"   - Calculation: {smoking_packs_per_day} packs/day * {smoking_years} years = {pack_years} pack-years.")
    print(" - Major risk factor: History of alcohol use disorder.")
    print("   - Both heavy smoking and alcohol are primary risk factors for Esophageal Squamous Cell Carcinoma (SCC).")
    print("-" * 50)

    print("Step 2: Analyzing Clinical Presentation and Test Results")
    print(" - Symptoms: Severe substernal chest pain and odynophagia (painful swallowing).")
    print(" - Labs: Elevated inflammatory markers (CRP, leukocytosis).")
    print(" - Imaging: Shows esophageal lumen narrowing and wall thickening.")
    print("   - This combination suggests a mass or an infiltrative process within the esophageal wall.")
    print(" - Endoscopy: Critically, the endoscopy is NORMAL, showing no erythema, ulcers, or plaques.")
    print("-" * 50)
    
    print("Step 3: Differential Diagnosis and Conclusion")
    print(" - The discrepancy between wall thickening on imaging and a normal-appearing mucosa on endoscopy is the key clue.")
    print(" - This pattern strongly suggests a submucosal or infiltrative lesion, where the disease grows within the wall without initially breaking the surface.")
    print(" - Infectious causes (A, E) and severe GERD (D) are effectively ruled out by the normal endoscopy, as they typically cause visible mucosal abnormalities (ulcers, plaques, erosions).")
    print(" - While adenocarcinoma (B) is a possibility, it is more associated with GERD/Barrett's esophagus. The risk factors of heavy smoking and alcohol use are far more specific to SCC.")
    print(" - Therefore, an infiltrative Squamous Cell Carcinoma (SCC) best explains all the findings: the risk factors, the symptoms, and the paradoxical combination of imaging and endoscopic results.")
    print("-" * 50)
    
    print("Most Likely Diagnosis: C. Esophageal squamous cell carcinoma")

diagnose_esophageal_condition()