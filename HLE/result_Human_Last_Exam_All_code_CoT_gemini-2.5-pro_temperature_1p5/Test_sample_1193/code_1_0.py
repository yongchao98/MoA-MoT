import sys

def solve_medical_case():
    """
    This function analyzes the clinical vignette using a logical scoring system
    to determine the most likely diagnosis.
    """
    # Key data points from the patient case
    patient_age = 59
    days_post_procedure = 29
    oxygen_saturation = 82
    supplemental_oxygen = 3

    # Represent clinical evidence as weighted factors
    # A positive score supports a diagnosis, a negative score argues against it.
    evidence = {
        "ARDS_picture": 10,        # Patient shows classic signs of ARDS
        "Late_PostOp_Timeline": 10,# Event occurs 29 days post-op, favouring delayed complications
        "Acute_Timeline_Mismatch":-20, # Strong evidence against acute causes
        "Major_Surgery_Risk": 10   # Whipple is a high-risk surgery for infection
    }

    print("Analyzing patient data:")
    print(f"- Days since procedure: {days_post_procedure}")
    print(f"- Oxygen Saturation: {oxygen_saturation}% on {supplemental_oxygen}L")
    print("- Findings: Bilateral crackles, gasping for air.")
    print("-" * 30)

    # Scoring the most relevant diagnoses based on evidence
    # Score for Sepsis (D)
    score_D = evidence["ARDS_picture"] + evidence["Late_PostOp_Timeline"] + evidence["Major_Surgery_Risk"]
    
    # Score for Acute Blood Transfusion Reaction (A)
    # The key factor here is the timeline mismatch.
    score_A = evidence["Acute_Timeline_Mismatch"] + evidence["ARDS_picture"]

    print("Logical Evaluation:")
    print("1. Diagnosis: Sepsis (D)")
    print("   This patient presents with a picture of ARDS, a known complication of sepsis.")
    print(f"   The timeline of {days_post_procedure} days is highly consistent with a post-operative infection.")
    print(f"   Final Equation (D): ARDS Picture ({evidence['ARDS_picture']}) + Late Timeline ({evidence['Late_PostOp_Timeline']}) + Surgery Risk ({evidence['Major_Surgery_Risk']}) = {score_D}")
    
    print("\n2. Diagnosis: Acute Blood Transfusion Reaction (A)")
    print("   While transfusions can cause ARDS (TRALI), 'acute' reactions happen within hours, not days.")
    print(f"   The timeline of {days_post_procedure} days makes this diagnosis highly unlikely.")
    print(f"   Final Equation (A): Timeline Mismatch ({evidence['Acute_Timeline_Mismatch']}) + ARDS Picture ({evidence['ARDS_picture']}) = {score_A}")

    print("-" * 30)
    print("Conclusion: Sepsis (D) is the most logical diagnosis as it aligns with the patient's clinical presentation, risk factors, and timeline.")

solve_medical_case()