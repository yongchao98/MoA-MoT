def solve_medical_case():
    """
    Analyzes a clinical case to find the most likely diagnosis.
    """
    # Step 1: Define patient's key clinical features
    patient_findings = {
        "History": "Faint cyanosis as a child, progressive dyspnea and fatigue in adulthood.",
        "Murmur Type": "Systolic Ejection",
        "Murmur Location": "Left Upper Sternal Border (LUSB)",
        "Murmur with Inspiration": "Increases",
        "Palpable Thrill": "Yes, at LUSB",
        "ECG Axis": "Left Axis Deviation (LAD)",
        "ECG Hypertrophy": "Right Ventricular Hypertrophy (RVH)"
    }

    # Step 2: Define features of possible diagnoses
    diagnoses = {
        "A. Ebstein anomaly": "Primarily causes tricuspid regurgitation (holosystolic murmur). ECG typically shows Right Bundle Branch Block and Right Axis Deviation.",
        "B. Patent ductus arteriosus": "Causes a continuous 'machine-like' murmur, not a systolic ejection murmur.",
        "C. Mitral valve prolapse": "Causes a mid-systolic click and/or late systolic murmur at the apex. It is a left-sided heart issue.",
        "D. Atrial septal defect": "Causes a systolic ejection 'flow' murmur at the LUSB due to increased volume across the pulmonary valve. Volume overload leads to RVH. An ostium primum type ASD is classically associated with Left Axis Deviation. A large, long-standing defect can present with symptoms in adulthood and a thrill.",
        "E. Hypertrophic cardiomyopathy with obstruction": "Murmur is at the left lower sternal border and DECREASES with inspiration. ",
        "F. Tricuspid stenosis": "Causes a diastolic murmur, not a systolic one.",
        "G. Ventricular septal defect": "Typically causes a holosystolic murmur at the left lower sternal border.",
        "H. None of these choices": "This is an option if no other choice fits."
    }

    # Step 3: Analyze the findings and evaluate choices
    print("Analyzing the Clinical Case:")
    print("-" * 30)
    print("Patient Presentation Summary:")
    for key, value in patient_findings.items():
        print(f"- {key}: {value}")
    print("-" * 30)

    print("Evaluating Answer Choices:")
    # The murmur characteristics (systolic ejection at LUSB, increases with inspiration) strongly point to a right-sided lesion, specifically increased flow or obstruction at the pulmonary valve.
    # The RVH is consistent with this.
    # The most unique and critical clue is the combination of RVH with Left Axis Deviation.
    
    print("\nKey Diagnostic Clues:")
    print("1. Murmur: The murmur is classic for pulmonic stenosis or a high-flow state across the pulmonary valve (like in an ASD).")
    print("2. Effect of Inspiration: The murmur increasing with inspiration (Carvallo's sign) confirms a right-sided heart pathology.")
    print("3. ECG Paradox: The combination of Right Ventricular Hypertrophy (from a right-sided problem) and Left Axis Deviation is highly specific.")
    print("\nAnalysis of Top Contender:")
    print("Choice D - Atrial septal defect (ASD):")
    print("- A flow murmur from an ASD matches the murmur's type and location.")
    print("- As a right-sided volume overload condition, it explains the RVH and the murmur's response to inspiration.")
    print("- An ostium primum ASD (a specific type of ASD) is a classic cause of Left Axis Deviation.")
    print("- A large, long-standing defect explains the adult-onset symptoms, faint childhood cyanosis (due to transient shunting), and a palpable thrill (from very high flow).")
    print("\nConclusion:")
    print("Atrial septal defect is the only choice that can account for the entire constellation of findings, especially the rare combination of RVH and LAD.")

    # Step 4: Final Answer
    final_answer = "D"
    print(f"\nThe most likely cause of this woman's murmur is D. Atrial septal defect.")

solve_medical_case()