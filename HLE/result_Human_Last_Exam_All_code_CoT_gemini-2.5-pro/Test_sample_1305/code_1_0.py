def solve_clinical_case():
    """
    This function analyzes the clinical vignette to determine the most likely diagnosis.
    """

    patient_findings = {
        "History": "Progressive shortness of breath, fatigue, childhood cyanosis",
        "Murmur Type": "Systolic ejection murmur",
        "Murmur Location": "Left upper sternal border (LUSB)",
        "Murmur Maneuver": "Increases with inspiration",
        "ECG": "Left Axis Deviation (LAD) and Right Ventricular Hypertrophy (RVH)"
    }

    print("Analyzing the patient's key findings:")
    for key, value in patient_findings.items():
        print(f"- {key}: {value}")

    print("\nEvaluating the answer choices:")
    print("A. Ebstein anomaly: Incorrect. Murmur is typically holosystolic at LLSB. ECG shows RBBB, not LAD.")
    print("B. Patent ductus arteriosus: Incorrect. Murmur is continuous ('machine-like').")
    print("C. Mitral valve prolapse: Incorrect. Left-sided lesion; murmur at apex.")
    print("D. Atrial septal defect: Correct. Let's analyze why:")
    print("  - A left-to-right shunt increases flow across the pulmonic valve, causing a systolic ejection murmur at the LUSB that increases with inspiration.")
    print("  - The volume overload on the right heart causes Right Ventricular Hypertrophy (RVH).")
    print("  - The specific ECG combination of Left Axis Deviation (LAD) with RVH is a classic sign of an ostium primum ASD, which is a type of atrial septal defect.")
    print("E. Hypertrophic cardiomyopathy: Incorrect. Murmur decreases with inspiration.")
    print("F. Tricuspid stenosis: Incorrect. Causes a diastolic murmur.")
    print("G. Ventricular septal defect: Incorrect. Murmur is typically holosystolic at the LLSB.")

    print("\nConclusion:")
    print("The combination of a systolic flow murmur at the LUSB and the specific ECG finding of LAD with RVH points strongly to an ostium primum Atrial Septal Defect.")
    
    final_answer_choice = "D"
    final_answer_text = "Atrial septal defect"
    print(f"The most likely cause is choice {final_answer_choice}: {final_answer_text}")

solve_clinical_case()