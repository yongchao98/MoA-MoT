def solve_diagnosis():
    """
    Analyzes the clinical case to determine the most likely diagnosis
    by scoring each option based on key findings.
    """
    # Define a dictionary for diagnoses and their likelihood scores
    diagnoses_scores = {
        'A. Streptococcal esophagitis': 0,
        'B. Esophageal adenocarcinoma': 0,
        'C. Esophageal squamous cell carcinoma': 0,
        'D. GERD': 0,
        'E. Herpes esophagitis': 0
    }

    # 1. Score based on major risk factors: 40+ pack-year smoking history and alcohol use disorder.
    # These are classic and very strong risk factors for Esophageal Squamous Cell Carcinoma (SCC).
    risk_factor_points = 50
    diagnoses_scores['C. Esophageal squamous cell carcinoma'] += risk_factor_points

    # 2. Score based on imaging vs. endoscopy findings.
    # Imaging shows wall thickening/narrowing, but endoscopy shows no mucosal lesions.
    # This pattern strongly suggests an infiltrative, submucosal tumor, characteristic of SCC.
    # It argues strongly against conditions that cause surface-level inflammation or ulcers.
    imaging_endoscopy_points_positive = 50
    imaging_endoscopy_points_negative = -50
    
    diagnoses_scores['C. Esophageal squamous cell carcinoma'] += imaging_endoscopy_points_positive
    diagnoses_scores['A. Streptococcal esophagitis'] += imaging_endoscopy_points_negative
    diagnoses_scores['D. GERD'] += imaging_endoscopy_points_negative
    diagnoses_scores['E. Herpes esophagitis'] += imaging_endoscopy_points_negative


    print("Calculating the likelihood score for the most probable diagnosis: Esophageal Squamous Cell Carcinoma (SCC).")
    print("The final score is the sum of points from key clinical findings.")
    print("\n--- Equation for SCC ---")
    
    initial_score = 0
    print(f"Initial Score: {initial_score}")
    print(f"Points from major risk factors (heavy smoking + alcohol): + {risk_factor_points}")
    print(f"Points from imaging (wall thickening) with negative endoscopy: + {imaging_endoscopy_points_positive}")
    
    final_scc_score = initial_score + risk_factor_points + imaging_endoscopy_points_positive
    
    print("---------------------------------")
    print(f"Total Score = {initial_score} + {risk_factor_points} + {imaging_endoscopy_points_positive} = {final_scc_score}")
    
    print("\nConclusion: Esophageal Squamous Cell Carcinoma is the diagnosis that best fits the combination of strong risk factors and the specific pattern of imaging and endoscopic findings.")

solve_diagnosis()