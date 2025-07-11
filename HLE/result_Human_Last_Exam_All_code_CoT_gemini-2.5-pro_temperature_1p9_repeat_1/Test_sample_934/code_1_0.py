def diagnose_esophageal_condition():
    """
    Analyzes a clinical vignette to determine the most likely diagnosis
    by scoring potential conditions based on the patient's data.
    """

    # --- Patient Data from Vignette ---
    # Strong risk factors for Esophageal Squamous Cell Carcinoma (SCC)
    heavy_smoking = True
    alcohol_use_disorder = True
    
    # Clinical Findings
    imaging_shows_wall_thickening = True
    endoscopy_is_negative_for_lesions = True # Key negative finding

    # --- Scoring System ---
    # We assign points based on how well the evidence supports each diagnosis.
    # Positive points for supporting findings, negative points for contradictory findings.

    # Score for Esophageal Squamous Cell Carcinoma (C)
    scc_score = 0
    # Major risk factors present
    scc_score_smoking = 5 if heavy_smoking else 0
    scc_score_alcohol = 5 if alcohol_use_disorder else 0
    # Imaging finding is classic for an infiltrative tumor
    scc_score_imaging = 3 if imaging_shows_wall_thickening else 0
    # Negative endoscopy is unusual but possible with submucosal/infiltrative SCC.
    # It strongly argues against other diagnoses that MUST have mucosal findings.
    # We penalize it slightly but recognize it doesn't rule out SCC.
    scc_score_endoscopy = -1 if endoscopy_is_negative_for_lesions else 0
    total_scc_score = scc_score_smoking + scc_score_alcohol + scc_score_imaging + scc_score_endoscopy

    # Scores for other conditions for comparison (showing why they are less likely)
    
    # B. Esophageal Adenocarcinoma: Smoking is a risk factor, but alcohol is not primary. Requires Barrett's/GERD history (absent).
    # Endoscopy would almost certainly show a mass. Negative endoscopy is strong evidence against it.
    adenocarcinoma_score = 2 + 3 - 5 # (smoking) + (imaging) - (negative endoscopy)
    
    # D. GERD: Does not cause significant wall thickening. Severe GERD causing this pain would show erosions on endoscopy.
    gerd_score = 1 - 3 - 5 # (risk factors) - (contradictory imaging) - (contradictory endoscopy)
    
    # A/E. Infectious Esophagitis (Strep/Herpes): Would show classic ulcers or plaques on endoscopy.
    # Negative endoscopy makes these extremely unlikely.
    infectious_score = 1 - 5 # (inflammation signs) - (contradictory endoscopy)
    
    
    # --- Conclusion ---
    scores = {
        "Esophageal squamous cell carcinoma": total_scc_score,
        "Esophageal adenocarcinoma": adenocarcinoma_score,
        "GERD": gerd_score,
        "Infectious esophagitis": infectious_score,
    }

    most_likely_diagnosis = max(scores, key=scores.get)
    
    print("Based on the patient's overwhelming risk factors and clinical findings, we can score the likelihood of each diagnosis.")
    print("The most likely diagnosis is Esophageal Squamous Cell Carcinoma (SCC).\n")
    print("Here is the scoring calculation for SCC:")
    
    # Print the "equation" as requested
    print(f"SCC Score = {scc_score_smoking} (for heavy smoking) + {scc_score_alcohol} (for heavy alcohol use) + {scc_score_imaging} (for wall thickening on imaging) + {scc_score_endoscopy} (for a non-contradictory endoscopy) = {total_scc_score}\n")

    print(f"Reasoning: The patient's history of heavy smoking and alcohol use are the strongest risk factors for SCC. The imaging findings of wall thickening and narrowing are classic for an infiltrative cancer. While the negative endoscopy is unusual, it can occur with submucosal tumors and argues more strongly against the other diagnoses, which would typically present with visible ulcers or plaques.")

diagnose_esophageal_condition()