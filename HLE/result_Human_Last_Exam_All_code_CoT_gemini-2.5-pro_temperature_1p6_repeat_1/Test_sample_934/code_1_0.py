import sys

def analyze_clinical_case():
    # --- Patient Data from Vignette ---
    age = 53
    pain_intensity = "10/10"
    cigarette_packs_per_day = 2
    smoking_years = 20

    # --- Key Findings ---
    risk_factors = ["Heavy alcohol use", "Heavy smoking"]
    imaging_findings = "Esophageal lumen narrowing and wall thickening"
    endoscopy_findings = "No signs of erythema, ulcers, plaques, or strictures"
    lab_findings = "Elevated C-reactive protein and leukocytosis (signs of inflammation)"

    # --- Diagnostic Reasoning ---
    print("Clinical Reasoning Steps:")
    print("-" * 25)

    # Step 1: Analyze major risk factors.
    # The two most significant risk factors for esophageal cancer are smoking and alcohol.
    # We can represent the smoking history numerically.
    pack_years = cigarette_packs_per_day * smoking_years
    print(f"1. The patient is a {age}-year-old with major risk factors for Esophageal Squamous Cell Carcinoma (SCC).")
    print(f"   - Smoking history calculation: {cigarette_packs_per_day} packs/day * {smoking_years} years = {pack_years} pack-years.")
    print(f"   - Significant alcohol use is also noted.")

    # Step 2: Evaluate the discrepancy between imaging and endoscopy.
    print("\n2. The key to the diagnosis is the conflict between imaging and endoscopy.")
    print(f"   - Imaging shows: '{imaging_findings}'. This points to a significant structural issue.")
    print(f"   - Endoscopy shows: '{endoscopy_findings}'. This rules out diseases that primarily affect the surface lining (mucosa).")

    # Step 3: Evaluate and eliminate other diagnoses.
    print("\n3. Eliminating other possibilities:")
    print("   - GERD, Herpes Esophagitis, and Streptococcal Esophagitis are unlikely because they typically cause visible mucosal abnormalities (redness, ulcers, plaques) which were absent.")
    print("   - Esophageal Adenocarcinoma is less likely, as its main risk factor is chronic GERD, not the smoking/alcohol profile seen here.")

    # Step 4: Conclude the most likely diagnosis.
    print("\n4. Conclusion:")
    print("   - Esophageal SCC can grow within the wall of the esophagus (submucosally) without initially breaking through the surface.")
    print("   - This pattern explains why imaging would show wall thickening and narrowing while the surface appears normal on endoscopy.")
    print(f"   - The patient's substernal chest pain, rated {pain_intensity}, is consistent with this invasive process.")

    print("-" * 25)
    print("Most Likely Diagnosis: C. Esophageal squamous cell carcinoma")

# Execute the analysis
analyze_clinical_case()