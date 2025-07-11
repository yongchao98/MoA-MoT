def analyze_ecg():
    """
    Analyzes the provided ECG findings and calculates the heart rate.
    """
    print("ECG Analysis:")
    print("1. Rhythm: The R-R intervals are variable, indicating an irregularly irregular rhythm. This is characteristic of Atrial Fibrillation.")
    print("2. QRS Complex: The QRS complexes are wide (> 0.12s) and their shape (morphology) changes from beat to beat.")
    print("3. Heart Rate: The rate is very fast (tachycardia).")

    # ECG paper has large squares representing 0.2 seconds.
    # Rate = 300 / (number of large squares between R waves)
    shortest_rr_interval_squares = 1.5
    longest_rr_interval_squares = 2.0

    rate_max = 300 / shortest_rr_interval_squares
    rate_min = 300 / longest_rr_interval_squares

    print(f"\nHeart Rate Calculation:")
    print(f"The shortest R-R interval is ~{shortest_rr_interval_squares} large squares.")
    print(f"  Calculation: 300 / {shortest_rr_interval_squares} = {rate_max:.0f} bpm")
    print(f"The longest R-R interval is ~{longest_rr_interval_squares} large squares.")
    print(f"  Calculation: 300 / {longest_rr_interval_squares} = {rate_min:.0f} bpm")
    print(f"The ventricular rate is dangerously fast, varying between {rate_min:.0f} and {rate_max:.0f} bpm.")

    print("\nConclusion:")
    print("The combination of an irregularly irregular rhythm (Atrial Fibrillation) with a very fast, wide, and morphologically varying QRS tachycardia is the classic presentation of Pre-excited Atrial Fibrillation (e.g., Wolff-Parkinson-White with AFib).")
    print("This occurs when atrial impulses bypass the AV node and travel down a faster accessory pathway.")

# Run the analysis
analyze_ecg()