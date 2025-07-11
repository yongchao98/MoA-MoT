def diagnose_heart_murmur():
    """
    Analyzes a clinical vignette to determine the most likely cause of a heart murmur.
    """

    # 1. Define the patient's key clinical findings from the vignette.
    patient_findings = {
        "Murmur Location": "Left upper sternal border",
        "Murmur Type": "Systolic ejection",
        "Effect of Inspiration": "Increases intensity (right-sided origin)",
        "EKG Axis": "Left axis deviation",
        "EKG Hypertrophy": "Right ventricular hypertrophy (RVH)",
        "Key Paradoxical Finding": "Left Axis Deviation + Right Ventricular Hypertrophy"
    }

    # 2. Explain the diagnostic reasoning step-by-step.
    print("Step 1: Analyzing the patient's key findings.")
    for key, value in patient_findings.items():
        print(f"- {key}: {value}")

    print("\nStep 2: Interpreting the findings.")
    print("The murmur's location (left upper sternal border) and type (systolic ejection) suggest increased blood flow across the pulmonic valve.")
    print("The murmur increasing with inspiration (Carvallo's sign) confirms its origin is on the right side of the heart.")
    print("The most critical clue is the EKG. Normally, RVH causes a Right Axis Deviation. The presence of a Left Axis Deviation with RVH is a paradoxical finding.")

    print("\nStep 3: Evaluating the diagnosis options.")
    print("This specific combination of a right-sided flow murmur and the paradoxical EKG (LAD + RVH) is a classic presentation for an ostium primum Atrial Septal Defect (ASD).")
    print("- The ASD causes a left-to-right shunt, leading to volume overload of the right ventricle, causing RVH and the pulmonic flow murmur.")
    print("- The primum type of ASD is an atrioventricular septal defect, which alters the heart's electrical conduction system, leading to a Left Axis Deviation.")

    print("\nStep 4: Conclusion.")
    print("Based on the unique combination of physical exam and EKG findings, Atrial Septal Defect is the most likely diagnosis among the choices.")

# Execute the analysis function
diagnose_heart_murmur()