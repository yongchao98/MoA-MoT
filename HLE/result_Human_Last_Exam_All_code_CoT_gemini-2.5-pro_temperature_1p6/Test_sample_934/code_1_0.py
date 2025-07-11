def diagnose_esophageal_condition():
    """
    Analyzes patient data to determine the most likely diagnosis.
    This function processes the clinical information from the vignette,
    evaluates risk factors and test results, and presents a reasoned conclusion.
    """

    # --- Patient Data from Vignette ---
    age_years = 53
    pain_intensity = 10
    cigarette_packs_per_day = 2
    smoking_duration_years = 20

    # Calculate pack-years, a key risk metric
    pack_years = cigarette_packs_per_day * smoking_duration_years

    print("--- Patient Clinical Summary ---")
    print(f"Patient Age: {age_years} years")
    print(f"Presenting Symptom: Substernal chest pain, intensity {pain_intensity}/10 with odynophagia.")
    print(f"Key Lab Findings: Elevated C-reactive protein and leukocytosis (signs of inflammation).")
    print("\n")

    print("--- Analysis of Risk Factors ---")
    print(f"Smoking History: {cigarette_packs_per_day} packs/day for >{smoking_duration_years} years, totaling over {pack_years} pack-years.")
    print("Alcohol History: Significant history of alcohol use disorder.")
    print("Conclusion on Risk: Heavy smoking and alcohol use are major risk factors for Esophageal Squamous Cell Carcinoma (SCC).")
    print("\n")

    print("--- Analysis of Diagnostic Findings ---")
    print("Imaging (CT/other): Shows esophageal lumen narrowing and wall thickening.")
    print("Endoscopy: Reveals NO visible signs of erythema, ulcers, or plaques (normal mucosa).")
    print("Diagnostic Inference: The discrepancy between abnormal imaging (thick wall) and normal endoscopy (clear surface) points to a submucosal infiltrative process.")
    print("\n")

    print("--- Differential Diagnosis Evaluation ---")
    print("A. Streptococcal/E. Herpes Esophagitis: Unlikely. These infections typically cause visible ulcers/plaques on endoscopy.")
    print("B. Esophageal Adenocarcinoma: Unlikely. Would almost certainly present as a visible mass/ulcer on endoscopy.")
    print("D. GERD: Unlikely. Does not typically cause significant wall thickening or the degree of inflammatory markers seen.")
    print("C. Esophageal Squamous Cell Carcinoma: Most Likely. This diagnosis aligns perfectly with the major risk factors (smoking, alcohol) and can present as an infiltrative, submucosal tumor, explaining the abnormal imaging despite a normal-appearing endoscopy.")

if __name__ == "__main__":
    diagnose_esophageal_condition()
<<<C>>>