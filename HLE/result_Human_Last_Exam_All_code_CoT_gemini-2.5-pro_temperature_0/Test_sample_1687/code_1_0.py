def diagnose_patient():
    """
    This script analyzes clinical findings from a case study to determine the most likely diagnosis.
    It uses a simple scoring system to weigh the significance of each finding for the differential diagnoses.
    """

    # --- Define weights for key clinical findings ---
    # These weights represent the clinical significance for this specific differential.
    luq_pain = 3
    left_shoulder_pain = 3  # Kehr's sign is highly suggestive of splenic injury
    rapid_hgb_drop = 2      # Indicates significant hemorrhage
    no_polypectomy = 0      # A factor to nullify diagnoses related to polypectomy

    # --- Calculate scores for each diagnosis based on the findings ---

    # A. Colonic perforation: Can cause pain and bleeding, but the specific pain location is less typical.
    score_perforation = 1 + rapid_hgb_drop

    # B. Lower GI bleeding: Explains blood loss but not the specific pain pattern.
    score_gi_bleed = 0 + rapid_hgb_drop

    # C. Splenic laceration: The classic triad of LUQ pain, shoulder pain, and hemorrhage after colonoscopy.
    score_splenic_laceration = luq_pain + left_shoulder_pain + rapid_hgb_drop

    # D. Postpolypectomy syndrome: Ruled out as no polypectomy was performed.
    score_postpolypectomy = (2 + 1) * no_polypectomy # The score will be 0

    diagnoses = {
        "A. Colonic perforation": score_perforation,
        "B. Lower GI bleeding": score_gi_bleed,
        "C. Splenic laceration": score_splenic_laceration,
        "D. Postpolypectomy syndrome": score_postpolypectomy
    }

    # Find the diagnosis with the highest score
    most_likely_diagnosis = max(diagnoses, key=diagnoses.get)
    highest_score = diagnoses[most_likely_diagnosis]

    print("Clinical Reasoning:")
    print("The patient's presentation with left upper quadrant pain, referred left shoulder pain (Kehr's sign), and hemodynamic shock after a difficult colonoscopy is a classic presentation of splenic injury.")
    print("\nScoring the diagnoses based on these findings:")
    
    # Print the equation for the highest-scoring diagnosis
    print(f"Most Likely Diagnosis Score (Splenic Laceration) = {luq_pain} (for LUQ pain) + {left_shoulder_pain} (for shoulder pain) + {rapid_hgb_drop} (for Hgb drop)")
    
    print(f"\nFinal Equation Result: {luq_pain} + {left_shoulder_pain} + {rapid_hgb_drop} = {highest_score}")

    print("\nConclusion:")
    print(f"The most likely diagnosis is '{most_likely_diagnosis}' with a score of {highest_score}.")

diagnose_patient()
<<<C>>>