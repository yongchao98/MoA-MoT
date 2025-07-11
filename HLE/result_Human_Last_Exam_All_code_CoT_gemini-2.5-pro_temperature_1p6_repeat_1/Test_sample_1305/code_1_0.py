def diagnose_murmur():
    """
    This function analyzes the clinical findings to determine the most likely diagnosis.
    """

    # Key Clinical Findings
    murmur_type = "Systolic ejection murmur"
    murmur_location = "Left upper sternal border"
    inspiration_effect = "Increases intensity"  # Carvallo's sign -> Right-sided lesion
    ecg_finding_1 = "Right ventricular hypertrophy (RVH)"
    ecg_finding_2 = "Left axis deviation (LAD)"

    # Analysis Steps
    print("Step 1: Analyze murmur characteristics.")
    print(f"A {murmur_type} at the {murmur_location} suggests increased flow across the pulmonic valve.")
    print(f"The murmur {inspiration_effect} with inspiration, pointing to a right-sided heart pathology.")
    print("-" * 30)

    print("Step 2: Analyze ECG findings.")
    print(f"The ECG shows {ecg_finding_1}, consistent with right-sided overload.")
    print(f"The ECG also shows {ecg_finding_2}.")
    print("The combination of RVH and LAD is a very specific finding, as RVH typically causes right axis deviation.")
    print("-" * 30)

    print("Step 3: Synthesize findings and conclude.")
    print("The constellation of findings (pulmonic flow murmur, RVH, and LAD) is a classic presentation for an ostium primum Atrial Septal Defect (ASD).")
    print("This defect causes a left-to-right shunt leading to the flow murmur and RVH, while the associated atrioventricular conduction defect causes the LAD.")
    print("-" * 30)

    # Final Answer
    final_diagnosis = "D. Atrial septal defect"
    print(f"The most likely cause is: {final_diagnosis}")

diagnose_murmur()