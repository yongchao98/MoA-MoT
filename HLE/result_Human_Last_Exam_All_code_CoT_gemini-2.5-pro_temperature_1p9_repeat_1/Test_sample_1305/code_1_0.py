def solve_clinical_vignette():
    """
    Analyzes the clinical vignette to determine the most likely diagnosis.
    """
    # Patient Presentation
    age = 31
    sex = "Female"
    symptoms = [
        "progressive shortness of breath on exertion",
        "mild fatigue",
        "occasional lightheadedness",
        "history of faint cyanosis as a child"
    ]

    # Physical Exam Findings
    murmur_type = "systolic ejection murmur"
    murmur_location = "left upper sternal border"
    murmur_dynamics = "increases in intensity with inspiration (Carvallo's sign)"

    # ECG Findings
    ecg_findings = ["left axis deviation (LAD)", "signs of right ventricular hypertrophy (RVH)"]

    print("Step 1: Analyzing the key findings.")
    print(f"  - Murmur: The murmur is a systolic ejection type at the left upper sternal border. This suggests increased blood flow across the pulmonic valve.")
    print(f"  - Dynamic Maneuver: The murmur increases with inspiration. This strongly points to a right-sided heart pathology, as inspiration increases venous return to the right heart.")
    print(f"  - ECG: The combination of Right Ventricular Hypertrophy (RVH) and Left Axis Deviation (LAD) is a crucial and paradoxical clue. Normally, RVH causes right axis deviation.")
    print("-" * 20)

    print("Step 2: Evaluating the answer choices.")
    # A. Ebstein anomaly: Typically causes tricuspid regurgitation (holosystolic murmur) and is not classically associated with LAD.
    # B. Patent ductus arteriosus: Causes a continuous 'machine-like' murmur, not a systolic ejection murmur.
    # C. Mitral valve prolapse: Left-sided event; murmur would decrease with inspiration.
    # E. HOCM: Left-sided event; murmur at a different location, decreases with inspiration.
    # F. Tricuspid stenosis: Causes a diastolic murmur, not systolic.
    # G. Ventricular septal defect: Holosystolic murmur, typically at the left LOWER sternal border. Does not cause the classic RVH + LAD combination.
    
    print("Step 3: Identifying the best fit.")
    print(f"  - The diagnosis must explain a flow murmur across the pulmonic valve (from a left-to-right shunt), signs of RVH (from right ventricular volume overload), and LAD.")
    print(f"  - An Atrial Septal Defect (ASD) causes a left-to-right shunt, leading to a flow murmur at the left upper sternal border and eventual RVH. This matches most of the signs.")
    print(f"  - The specific ECG finding of RVH with LAD is pathognomonic for a specific type of ASD: an ostium primum ASD (a form of endocardial cushion defect). This type of defect is associated with a conduction abnormality that causes LAD.")
    print(f"  - Since ostium primum ASD is a subtype of Atrial Septal Defect, choice D is the most fitting diagnosis.")
    print("-" * 20)

    final_answer_choice = "D"
    final_answer_text = "Atrial septal defect"

    print(f"Conclusion: The constellation of findings, especially the paradoxical ECG, is classic for an ostium primum Atrial Septal Defect. Therefore, 'Atrial septal defect' is the most likely cause.")
    
# Execute the function to print the reasoning.
solve_clinical_vignette()