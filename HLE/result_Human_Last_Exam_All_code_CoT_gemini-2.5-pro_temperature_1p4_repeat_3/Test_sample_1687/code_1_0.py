def diagnose_colonoscopy_complication():
    """
    This script analyzes clinical findings from a post-colonoscopy case
    to determine the most likely diagnosis. It uses a scoring system
    to weigh the evidence for each possibility.
    """

    # --- Patient's Key Clinical Findings ---
    # We represent the presence of a key finding with a value of 1.
    patient_findings = {
        "Left Upper Quadrant (LUQ) Pain": 1,
        "Left Shoulder Pain (Kehr's Sign)": 1,
        "Severe Hemorrhage / Shock": 1,
        "Fact: No Polypectomy Performed": 1,
    }

    # --- Diagnostic Scoring Matrix ---
    # Each diagnosis is given points based on how strongly a finding supports it.
    # High score (2): Classic sign for the diagnosis.
    # Medium score (1): Consistent with the diagnosis.
    # Neutral score (0): Not typically associated.
    # Negative score (-10): Contradicts or rules out the diagnosis.
    diagnostic_criteria = {
        "A. Colonic perforation": {
            "Left Upper Quadrant (LUQ) Pain": 1,
            "Left Shoulder Pain (Kehr's Sign)": 0,
            "Severe Hemorrhage / Shock": 1, # Possible but more often septic shock
            "Fact: No Polypectomy Performed": 0,
        },
        "B. Lower GI bleeding": {
            "Left Upper Quadrant (LUQ) Pain": 0, # Pain is not a primary feature
            "Left Shoulder Pain (Kehr's Sign)": 0,
            "Severe Hemorrhage / Shock": 2, # Is a hemorrhagic condition
            "Fact: No Polypectomy Performed": 0,
        },
        "C. Splenic laceration": {
            "Left Upper Quadrant (LUQ) Pain": 2, # Classic location of spleen
            "Left Shoulder Pain (Kehr's Sign)": 2, # Classic referred pain from diaphragm irritation
            "Severe Hemorrhage / Shock": 2, # Spleen is highly vascular, causing massive bleed
            "Fact: No Polypectomy Performed": 0,
        },
        "D. Postpolypectomy syndrome": {
            "Left Upper Quadrant (LUQ) Pain": 1,
            "Left Shoulder Pain (Kehr's Sign)": 0,
            "Severe Hemorrhage / Shock": 0, # Not a hemorrhagic process
            "Fact: No Polypectomy Performed": -10, # Ruled out if no polypectomy was done
        },
    }

    results = {}
    best_diagnosis = ""
    highest_score = -99

    print("--- Diagnostic Analysis ---\n")

    for diagnosis, criteria in diagnostic_criteria.items():
        score = 0
        equation_str = f"Calculating score for {diagnosis}:\n"
        
        # Calculate score by multiplying patient finding (0 or 1) with criteria weight
        luq_pain_score = patient_findings["Left Upper Quadrant (LUQ) Pain"] * criteria["Left Upper Quadrant (LUQ) Pain"]
        shoulder_pain_score = patient_findings["Left Shoulder Pain (Kehr's Sign)"] * criteria["Left Shoulder Pain (Kehr's Sign)"]
        hemorrhage_score = patient_findings["Severe Hemorrhage / Shock"] * criteria["Severe Hemorrhage / Shock"]
        polypectomy_score = patient_findings["Fact: No Polypectomy Performed"] * criteria["Fact: No Polypectomy Performed"]
        
        score = luq_pain_score + shoulder_pain_score + hemorrhage_score + polypectomy_score
        
        equation_str += (f"  LUQ Pain Contribution: {luq_pain_score}\n"
                         f"  Left Shoulder Pain Contribution: {shoulder_pain_score}\n"
                         f"  Hemorrhage/Shock Contribution: {hemorrhage_score}\n"
                         f"  No Polypectomy Factor: {polypectomy_score}\n")

        results[diagnosis] = score
        equation_str += f"  Total Score = {luq_pain_score} + {shoulder_pain_score} + {hemorrhage_score} + {polypectomy_score} = {score}\n"
        print(equation_str)

        if score > highest_score:
            highest_score = score
            best_diagnosis = diagnosis

    print("--- Conclusion ---")
    print(f"The diagnosis with the highest likelihood score is: {best_diagnosis} (Score: {highest_score})")
    print("\nReasoning: Splenic laceration uniquely explains the combination of Left Upper Quadrant pain, referred left shoulder pain (Kehr's sign), and profound hypovolemic shock from internal hemorrhage following a difficult colonoscopy.")
    
# Execute the diagnostic analysis
diagnose_colonoscopy_complication()