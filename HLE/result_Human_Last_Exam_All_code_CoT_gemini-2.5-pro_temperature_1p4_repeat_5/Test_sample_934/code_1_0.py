def solve_medical_case():
    """
    Analyzes the clinical vignette to determine the most likely diagnosis.
    This function will print the step-by-step reasoning based on the provided patient data.
    """

    # --- Step 1: Quantify Key Risk Factors ---
    print("--- Step 1: Analyzing Patient Data and Risk Factors ---")
    packs_per_day = 2
    years_smoking = 20
    # The "equation" to calculate pack-years:
    pack_years = packs_per_day * years_smoking
    
    print(f"The patient is a 53-year-old woman.")
    print(f"She has a significant smoking history of {packs_per_day} packs/day for {years_smoking} years.")
    print(f"Calculating smoking pack-years: {packs_per_day} * {years_smoking} = {pack_years} pack-years.")
    print("Other major risk factors include a history of alcohol use disorder.")
    print("Conclusion for Step 1: The patient's profile (age, heavy smoking, alcohol use) presents classic and strong risk factors for Esophageal Squamous Cell Carcinoma (SCC).\n")

    # --- Step 2: Evaluating Clinical and Diagnostic Findings ---
    print("--- Step 2: Evaluating Clinical and Diagnostic Findings ---")
    print("Symptoms: Substernal chest pain (10/10) and odynophagia (pain with swallowing) are consistent with esophageal cancer.")
    print("Imaging: 'Esophageal lumen narrowing wall thickening' is highly suggestive of an infiltrative tumor.")
    print("Endoscopy: This is a crucial finding. The test was negative for erythema, ulcers, plaques, or strictures.")
    print("Conclusion for Step 2: The absence of mucosal findings on endoscopy makes infectious esophagitis (Herpes, Strep) or severe GERD very unlikely. An infiltrative SCC growing within the esophageal wall can explain the imaging findings while the surface mucosa appears normal.\n")

    # --- Step 3: Final Diagnosis ---
    print("--- Step 3: Reaching a Final Diagnosis ---")
    print("The diagnosis must explain the combination of overwhelming risk factors, symptoms, and the specific pattern of diagnostic findings (positive imaging, negative endoscopy).")
    print("Esophageal Squamous Cell Carcinoma is the only option that fits all these criteria perfectly.\n")
    
    print("Most Likely Diagnosis: Esophageal squamous cell carcinoma")

# Execute the analysis
solve_medical_case()