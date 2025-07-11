def diagnose_ecg():
    """
    This function analyzes ECG characteristics to determine the diagnosis.
    """
    print("Starting ECG analysis...")
    
    # Step 1: Analyze the Rhythm (R-R interval)
    print("\n--- Step 1: Rhythm Analysis ---")
    rhythm = "Irregularly Irregular"
    print(f"Observation: The R-R intervals are not constant.")
    print(f"Conclusion: The rhythm is {rhythm}.")
    
    # Step 2: Analyze the Heart Rate
    print("\n--- Step 2: Heart Rate Analysis ---")
    shortest_rr_interval_squares = 1.5
    rate = 300 / shortest_rr_interval_squares
    print(f"Observation: The shortest R-R interval is approximately {shortest_rr_interval_squares} large squares.")
    print(f"Calculation: Rate = 300 / {shortest_rr_interval_squares} = {rate} bpm.")
    print("Conclusion: This is a very rapid Tachycardia (>150 bpm).")

    # Step 3: Analyze the QRS Width
    print("\n--- Step 3: QRS Width Analysis ---")
    qrs_width_seconds = 0.16 # estimated > 0.12s
    print(f"Observation: The QRS duration is > 0.12 seconds (3 small squares).")
    print("Conclusion: This is a Wide-Complex Tachycardia.")

    # Step 4: Synthesize Findings and reach a conclusion
    print("\n--- Step 4: Synthesis and Final Diagnosis ---")
    print(f"Summary of findings: A very rapid ({int(rate)} bpm), {rhythm}, Wide-Complex Tachycardia.")
    print("Additionally, the QRS morphology shows beat-to-beat variation.")
    print("\nThis specific combination of features is classic for atrial fibrillation impulses conducting rapidly down an accessory pathway.")
    print("This bypasses the AV node's normal rate-limiting function, leading to the extremely high ventricular rate.")
    
    final_diagnosis = "Pre-excited Atrial Fibrillation"
    print(f"\nFinal Diagnosis: {final_diagnosis}")

diagnose_ecg()