def solve_clinical_case():
    """
    Analyzes a clinical vignette to determine the most likely diagnosis.
    """
    # 1. Key findings from the patient's case
    patient_findings = {
        "History": "31-year-old woman, history of faint cyanosis as a child",
        "Symptoms": "Progressive shortness of breath, fatigue, lightheadedness",
        "Murmur_Type": "Systolic ejection murmur",
        "Murmur_Location": "Left upper sternal border (pulmonic area)",
        "Murmur_Maneuver": "Increases in intensity with inspiration (Carvallo's sign)",
        "ECG_Findings": ["Left axis deviation (LAD)", "Right ventricular hypertrophy (RVH)"]
    }

    # 2. Evaluation of Findings
    print("Thinking Process:")
    print("1. The murmur is a 'systolic ejection' type at the 'left upper sternal border'. This location corresponds to the pulmonic valve area. This suggests either pulmonic stenosis or increased blood flow across the pulmonic valve (relative stenosis).")
    print("2. The murmur 'increases with inspiration'. This is known as Carvallo's sign and strongly indicates a right-sided heart pathology, as inspiration increases venous return to the right heart.")
    print("3. The ECG shows a combination of 'Right Ventricular Hypertrophy (RVH)' and 'Left Axis Deviation (LAD)'.")
    print("   - RVH is consistent with the right-sided volume/pressure overload suggested by the murmur.")
    print("   - LAD is the key distinguishing feature. Most causes of RVH (like typical ASD or pulmonic stenosis) lead to Right Axis Deviation. The combination of RVH and LAD is very specific.")
    print("\nEvaluating the choices:")

    # 3. Comparing with answer choices
    analysis = {
        "A. Ebstein anomaly": "Typically causes a holosystolic murmur of tricuspid regurgitation and ECG shows Right Axis Deviation. Mismatch.",
        "B. Patent ductus arteriosus": "Causes a continuous 'machine-like' murmur, not a systolic ejection murmur. Mismatch.",
        "C. Mitral valve prolapse": "Causes a mid-systolic click and late systolic murmur at the apex. Mismatch.",
        "D. Atrial septal defect": "Causes increased flow across the pulmonic valve, leading to a systolic ejection murmur at the left upper sternal border and RVH. A specific type, ostium primum ASD, is classically associated with the paradoxical combination of Left Axis Deviation and RVH. This is a strong match.",
        "E. Hypertrophic cardiomyopathy": "Murmur is typically at the left lower sternal border and DECREASES with inspiration. Mismatch.",
        "F. Tricuspid stenosis": "Causes a diastolic murmur, not systolic. Mismatch.",
        "G. Ventricular septal defect": "Typically causes a holosystolic murmur at the left lower sternal border. An inlet VSD can cause LAD, but the classic ASD picture fits the murmur better. Less likely.",
    }

    for choice, reason in analysis.items():
        print(f"- {choice}: {reason}")

    # 4. Conclusion
    print("\nConclusion:")
    print("The combination of a systolic ejection murmur at the pulmonic area, signs of right heart volume overload (murmur increases with inspiration, RVH), and a paradoxical Left Axis Deviation on ECG is highly characteristic of an ostium primum Atrial Septal Defect (ASD). Therefore, Atrial Septal Defect is the most likely diagnosis among the choices.")

solve_clinical_case()