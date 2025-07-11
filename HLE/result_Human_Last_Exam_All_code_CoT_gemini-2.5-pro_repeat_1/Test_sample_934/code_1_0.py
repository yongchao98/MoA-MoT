def diagnose_esophageal_condition():
    """
    Analyzes a clinical case to determine the most likely diagnosis
    and prints the reasoning.
    """
    print("Analyzing the clinical case based on the provided information.")
    print("------------------------------------------------------------\n")

    # Step 1: Deconstruct the Patient Presentation
    print("Step 1: Deconstructing the Patient Presentation")
    print("- Patient: 53-year-old woman.")
    print("- Symptoms: 10/10 substernal chest pain and odynophagia (painful swallowing).")
    print("- Key Risk Factors: Heavy smoking (2 packs/day for 20 years) and alcohol use disorder. These are the two most significant risk factors for Esophageal Squamous Cell Carcinoma (SCC).")
    print("- Imaging Findings: Esophageal lumen narrowing and wall thickening. This suggests an infiltrative process, such as a tumor growing within the wall of the esophagus.")
    print("- Endoscopic Findings: The mucosa is normal, with no erythema, ulcers, plaques, or strictures. This is a critical finding.\n")

    # Step 2: Evaluating the Differential Diagnosis
    print("Step 2: Evaluating the Differential Diagnosis")
    print("- A. Streptococcal esophagitis, D. GERD, and E. Herpes esophagitis are all unlikely. These conditions cause inflammation and damage to the esophageal lining (mucosa). A normal endoscopic exam effectively rules them out as the primary cause of such severe symptoms and imaging findings.")
    print("- B. Esophageal adenocarcinoma is less likely. Its main risk factor is chronic GERD leading to Barrett's esophagus. The patient's risk profile points away from this diagnosis.")
    print("- C. Esophageal squamous cell carcinoma (SCC) is the most likely diagnosis. It aligns perfectly with all the data:")
    print("  1. The patient's risk factors (smoking, alcohol) are classic for SCC.")
    print("  2. The symptoms (pain, odynophagia) and imaging findings (wall thickening) are consistent with an invasive cancer.")
    print("  3. Most importantly, SCC can begin and grow within the submucosa (the layer beneath the lining). This causes the wall to thicken and the lumen to narrow without any visible changes on the surface, explaining the normal endoscopy.\n")

    # Step 3: Conclusion
    print("Step 3: Conclusion")
    print("The combination of classic risk factors, suggestive symptoms and imaging, and a normal endoscopy points strongly to an infiltrative submucosal tumor. Esophageal Squamous Cell Carcinoma is the most fitting diagnosis.")

diagnose_esophageal_condition()