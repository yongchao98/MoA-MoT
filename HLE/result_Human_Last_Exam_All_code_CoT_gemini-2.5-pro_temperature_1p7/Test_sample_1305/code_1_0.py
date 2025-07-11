import sys

def solve_clinical_case():
    """
    This function analyzes a clinical vignette to determine the most likely diagnosis.
    It prints the step-by-step reasoning.
    """
    
    # Step 1: Deconstruct the clinical vignette.
    print("Step 1: Deconstruct the clinical vignette.")
    print("Patient: 31-year-old woman.")
    print("Symptoms: Progressive shortness of breath on exertion, fatigue, lightheadedness, history of faint childhood cyanosis.")
    print("Physical Exam: Systolic ejection murmur at the left upper sternal border (LUSB), radiating to the left clavicle. Murmur increases in intensity with inspiration (Carvallo's sign). Palpable thrill.")
    print("ECG: Left axis deviation (LAD) and signs of right ventricular hypertrophy (RVH).\n")

    # Step 2: Analyze the key findings.
    print("Step 2: Analyze the key findings.")
    print("- The murmur characteristics (systolic ejection at LUSB, increasing with inspiration) strongly suggest pathology involving the right side of the heart, specifically increased blood flow across the pulmonic valve.")
    print("- The ECG combination of Left Axis Deviation (LAD) and Right Ventricular Hypertrophy (RVH) is a paradoxical but highly specific finding. Normally, RVH leads to Right Axis Deviation. This specific combination points strongly towards a particular type of congenital heart defect.\n")

    # Step 3: Evaluate the answer choices based on the key findings.
    print("Step 3: Evaluate the answer choices.")
    print("A. Ebstein anomaly: Incorrect. Typically causes a holosystolic murmur of tricuspid regurgitation and is not associated with the LAD + RVH combination.")
    print("B. Patent ductus arteriosus: Incorrect. Causes a continuous 'machinery' murmur.")
    print("C. Mitral valve prolapse: Incorrect. This is a left-sided valve issue, and its murmur (mid-systolic click) does not increase with inspiration.")
    print("D. Atrial septal defect (ASD): Correct. A large ASD causes a left-to-right shunt. This leads to volume overload of the right ventricle (causing RVH) and increased flow across the pulmonic valve (causing the systolic ejection 'flow murmur'). While a secundum ASD (the most common type) causes Right Axis Deviation, an ostium primum ASD (a type of atrioventricular septal defect) is classically associated with the specific ECG finding of LAD plus RVH. Therefore, the general category of ASD is the best fit for the complete clinical picture.")
    print("E. Hypertrophic cardiomyopathy: Incorrect. The murmur of HOCM typically decreases with maneuvers that increase venous return, like inspiration.")
    print("F. Tricuspid stenosis: Incorrect. This would cause a diastolic murmur.")
    print("G. Ventricular septal defect: Incorrect. Typically causes a holosystolic murmur at the lower left sternal border.\n")

    # Step 4: Conclusion.
    print("Step 4: Conclusion.")
    print("The constellation of a systolic flow murmur over the pulmonic area, evidence of right ventricular volume overload (RVH), and the specific ECG pattern of left axis deviation is classic for an ostium primum atrial septal defect. Therefore, 'Atrial septal defect' is the most likely cause of this woman's murmur.")

# Execute the analysis
solve_clinical_case()