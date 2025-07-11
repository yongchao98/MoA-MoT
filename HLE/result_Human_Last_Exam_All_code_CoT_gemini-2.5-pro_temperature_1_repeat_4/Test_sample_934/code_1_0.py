def diagnose_patient():
    """
    Analyzes the patient's clinical data to determine the most likely diagnosis
    by weighing the evidence for each potential condition.
    """

    # Patient Data
    patient_age = 53
    risk_factors = ["Heavy Smoking (>20 years)", "Alcohol Use Disorder"]
    symptoms = ["10/10 Substernal Chest Pain", "Pain with Swallowing (Odynophagia)"]
    lab_findings = ["Elevated C-reactive protein", "Leukocytosis"]
    imaging_findings = ["Esophageal lumen narrowing", "Wall thickening"]
    endoscopy_finding = "No visible erythema, ulcers, or plaques"

    print("Analyzing the diagnostic evidence:")
    print("=" * 50)
    print(f"Patient Profile: A {patient_age}-year-old woman.")
    print(f"Key Risk Factors: {risk_factors[0]} + {risk_factors[1]}.")
    print("These are major risk factors for Esophageal Squamous Cell Carcinoma (SCC).\n")

    print("Evaluating symptoms and findings against possible diagnoses:\n")

    # Analysis of GERD and Infectious Esophagitis
    print("1. GERD / Herpes Esophagitis / Streptococcal Esophagitis:")
    print(f"   - The key finding against these is: '{endoscopy_finding}'.")
    print("   - These conditions typically cause visible mucosal changes (ulcers, redness, plaques). Their absence makes these diagnoses unlikely.\n")

    # Analysis of Esophageal Adenocarcinoma
    print("2. Esophageal Adenocarcinoma:")
    print("   - This is less likely because the patient's primary risk factors (Smoking + Alcohol) point more strongly towards SCC.")
    print("   - This cancer also typically presents as a visible mass or ulcer on endoscopy.\n")

    # Analysis of Esophageal Squamous Cell Carcinoma (SCC)
    print("3. Esophageal Squamous Cell Carcinoma (SCC):")
    print("   - This diagnosis aligns with all the evidence.")
    print(f"   - Patient has the classic risk factors: {risk_factors[0]} and {risk_factors[1]}.")
    print(f"   - Imaging findings ({', '.join(imaging_findings)}) are explained by a tumor growing in the wall.")
    print(f"   - The normal endoscopy ('{endoscopy_finding}') is explained by a submucosal or infiltrative tumor that has not broken through the surface lining.\n")

    print("=" * 50)
    print("Final Diagnostic Conclusion:")
    print("The combination of:")
    print("Risk Factors (Smoking + Alcohol)")
    print("+ Imaging (Wall Thickening)")
    print("+ Normal Endoscopy (ruling out mucosal disease)")
    print("-----------------------------------------------------")
    print("= Most Likely Diagnosis: C. Esophageal Squamous Cell Carcinoma")

diagnose_patient()