def solve_clinical_case():
    """
    Analyzes a clinical vignette to determine the most likely diagnosis.
    """
    print("Step 1: Deconstructing the clinical vignette.")
    print("  - Patient: 31-year-old woman.")
    print("  - Symptoms: Shortness of breath, fatigue, lightheadedness (suggesting a significant cardiac issue).")
    print("  - History: Faint cyanosis as a child (suggests a congenital heart defect).")
    print("  - Physical Exam (Murmur): Systolic ejection murmur at the left upper sternal border (LUSB), increases with inspiration (Carvallo's sign).")
    print("  - Electrocardiogram (ECG): Left axis deviation (LAD) and signs of right ventricular hypertrophy (RVH).")
    print("\n")

    print("Step 2: Analyzing the key findings.")
    print("  - Murmur Analysis: A systolic ejection murmur at the LUSB is characteristic of increased blood flow across the pulmonary valve. The murmur increasing with inspiration (a right-sided event) confirms a pathology involving the right side of the heart.")
    print("  - ECG Analysis: The combination of Left Axis Deviation (LAD) with Right Ventricular Hypertrophy (RVH) is a very specific and relatively uncommon finding. RVH is expected with right-sided pressure or volume overload. However, these conditions (like pulmonic stenosis or a typical secundum ASD) usually cause Right Axis Deviation, not Left.")
    print("\n")

    print("Step 3: Evaluating the answer choices based on the key findings.")
    print("  - A. Ebstein anomaly: Causes tricuspid regurgitation (holosystolic murmur) and typically right axis deviation.")
    print("  - B. Patent ductus arteriosus: Causes a continuous 'machine-like' murmur, not a systolic ejection murmur.")
    print("  - C. Mitral valve prolapse: A left-sided issue with a murmur (mid-systolic click) best heard at the apex.")
    print("  - D. Atrial septal defect (ASD): This is a strong candidate. An ASD causes a left-to-right shunt, leading to volume overload of the right ventricle and increased flow across the pulmonary valve, which explains the murmur and the RVH. Crucially, a specific type of ASD, the 'ostium primum' ASD (which is a form of atrioventricular septal defect), is classically associated with Left Axis Deviation due to abnormal cardiac conduction pathways. This diagnosis perfectly fits the entire clinical picture.")
    print("  - E. Hypertrophic cardiomyopathy: Murmur typically increases with maneuvers that decrease venous return (like Valsalva), the opposite of what is seen here.")
    print("  - F. Tricuspid stenosis: Causes a diastolic murmur, not systolic.")
    print("  - G. Ventricular septal defect: Typically causes a holosystolic murmur at the lower left sternal border. An isolated VSD does not typically cause LAD.")
    print("\n")
    
    print("Step 4: Conclusion.")
    print("The combination of a right-sided volume overload pattern (systolic ejection murmur at LUSB from pulmonary flow, RVH) with Left Axis Deviation on ECG is a classic presentation for an ostium primum Atrial Septal Defect.")
    print("Therefore, the most likely cause is an Atrial Septal Defect.")
    print("\n")
    
    final_answer = "D"
    print(f"Final Answer: {final_answer}")

solve_clinical_case()