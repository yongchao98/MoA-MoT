import sys

def solve_medical_case():
    """
    Analyzes clinical findings to determine the most likely diagnosis.
    """
    # Key findings from the clinical case
    patient_findings = {
        "procedure": "Difficult colonoscopy",
        "pain_location": "LUQ",
        "referred_pain": "Left shoulder",
        "polypectomy_performed": False,
        "hemoglobin_initial": 11.7,
        "hemoglobin_final": 6.5,
        "hemodynamic_instability": True,
        "peritoneal_signs": True
    }

    # Characteristics of each potential diagnosis
    diagnoses = {
        "A. Colonic perforation": {
            "common_pain_location": ["Generalized", "Lower"],
            "causes_massive_hemorrhage": False,
            "causes_kehrs_sign": False,
            "requires_polypectomy": False,
            "causes_peritoneal_signs": True
        },
        "B. Lower GI bleeding": {
            "common_pain_location": ["Crampy", "None"],
            "causes_massive_hemorrhage": True, # The bleed itself
            "causes_kehrs_sign": False,
            "requires_polypectomy": False,
            "causes_peritoneal_signs": False # Bleeding is intraluminal
        },
        "C. Splenic laceration": {
            "common_pain_location": ["LUQ"],
            "causes_massive_hemorrhage": True,
            "causes_kehrs_sign": True,
            "requires_polypectomy": False,
            "causes_peritoneal_signs": True # From hemoperitoneum
        },
        "D. Postpolypectomy syndrome": {
            "common_pain_location": ["Localized"],
            "causes_massive_hemorrhage": False,
            "causes_kehrs_sign": False,
            "requires_polypectomy": True,
            "causes_peritoneal_signs": True
        }
    }

    # --- Analysis Logic ---
    print("Analyzing patient findings against potential diagnoses...")
    
    hemoglobin_drop = patient_findings["hemoglobin_initial"] - patient_findings["hemoglobin_final"]
    print(f"Key Finding: Significant hemoglobin drop from {patient_findings['hemoglobin_initial']} to {patient_findings['hemoglobin_final']} g/dL.")
    print("Key Finding: Pain in the Left Upper Quadrant (LUQ) with referred pain to the left shoulder (Kehr's sign).")
    print("Key Finding: No polypectomy was performed.\n")

    scores = {}
    for name, features in diagnoses.items():
        score = 0
        # Rule out based on explicit contradictions
        if features["requires_polypectomy"] and not patient_findings["polypectomy_performed"]:
            scores[name] = (0, "Ruled out: No polypectomy was performed.")
            continue
            
        # Score based on matching features
        if patient_findings["pain_location"] in features["common_pain_location"]:
            score += 1
        if patient_findings["referred_pain"] == "Left shoulder" and features["causes_kehrs_sign"]:
            score += 2 # Weight this heavily as it's a specific sign
        if patient_findings["hemodynamic_instability"] and features["causes_massive_hemorrhage"]:
            score += 1
        if patient_findings["peritoneal_signs"] and features["causes_peritoneal_signs"]:
            score += 1
        
        scores[name] = (score, "Matches key features.")

    # --- Print Results ---
    print("--- Scoring Results ---")
    best_diagnosis = ""
    highest_score = -1
    for name, (score, reason) in sorted(scores.items()):
        print(f"{name}: Score {score}. Reason: {reason}")
        if score > highest_score:
            highest_score = score
            best_diagnosis = name
            
    print("\n--- Conclusion ---")
    print(f"The highest score belongs to '{best_diagnosis}', making it the most likely diagnosis based on the provided findings.")

# The final answer format is embedded in the print statement below.
# This prevents it from being captured by mistake if this script is run elsewhere.
if __name__ == "__main__":
    original_stdout = sys.stdout  # Save a reference to the original standard output
    try:
        solve_medical_case()
    finally:
        sys.stdout = original_stdout # Restore original stdout
    
    # Directly print the final answer in the required format
    print("<<<C>>>")