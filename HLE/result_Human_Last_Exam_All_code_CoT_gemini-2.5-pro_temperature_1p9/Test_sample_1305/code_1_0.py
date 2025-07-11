import textwrap

def solve_medical_case():
    """
    Analyzes a clinical vignette to determine the most likely diagnosis.
    """

    explanation = [
        "1. Analyzing the Patient's Presentation:",
        "   - History: 31-year-old woman with lifelong symptoms (childhood cyanosis) that are now progressing (shortness of breath, fatigue). This suggests a congenital heart defect.",
        "   - Murmur: A systolic ejection murmur at the left upper sternal border that increases with inspiration. The location points to the pulmonary valve area. The increase with inspiration (Carvallo's sign) indicates a right-sided heart event, as inspiration increases blood flow to the right side of the heart.",
        "   - ECG: Shows signs of Right Ventricular Hypertrophy (RVH) and Left Axis Deviation (LAD). RVH is expected with a right-sided lesion causing volume or pressure overload. However, the combination of RVH with LAD is a very specific finding.",
        "\n2. Evaluating the Options:",
        "   - A, F (Ebstein, Tricuspid Stenosis): These are right-sided lesions, but Tricuspid Stenosis causes a diastolic murmur, and Ebstein anomaly classically causes a holosystolic murmur of tricuspid regurgitation, not a systolic ejection murmur.",
        "   - B (Patent Ductus Arteriosus): Causes a continuous 'machine-like' murmur, not a systolic ejection murmur.",
        "   - C, E (Mitral Prolapse, HOCM): These are left-sided lesions, and their murmurs would typically decrease, not increase, with inspiration.",
        "   - G (Ventricular Septal Defect): Classically causes a harsh, holosystolic murmur, typically at the left lower sternal border.",
        "\n3. Conclusion:",
        "   - D (Atrial Septal Defect): This diagnosis fits all the findings perfectly. An ASD is a left-to-right shunt that causes significant volume overload on the right heart, leading to RVH. This large volume of blood flowing across the pulmonary valve creates a 'relative pulmonic stenosis', resulting in a systolic ejection murmur exactly where described (left upper sternal border). Importantly, the specific subtype known as an 'ostium primum' ASD is characteristically associated with the unique ECG finding of Left Axis Deviation coexisting with Right Ventricular Hypertrophy.",
        "   - Therefore, an Atrial Septal Defect (specifically, ostium primum type) is the most likely cause of this woman's murmur and overall clinical picture."
    ]

    for line in explanation:
        # Use textwrap to format the lines for better readability
        print(textwrap.fill(line, width=100))
    
    # Final answer is D
    final_answer = "D"
    print(f"\n<<<D>>>")

solve_medical_case()