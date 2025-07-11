def analyze_ecg():
    """
    Analyzes the provided ECG and explains the diagnosis step-by-step.
    """
    print("ECG Analysis:")
    print("Step 1: Rhythm Analysis")
    print("Observation: The R-R intervals are not constant. They vary significantly from one beat to the next.")
    print("Conclusion: The rhythm is irregularly irregular, a hallmark of Atrial Fibrillation (AFib).\n")

    print("Step 2: Heart Rate Calculation")
    # In the ECG, the shortest R-R interval is about 7 small squares.
    # Standard ECG paper speed: 1 small square = 0.04 seconds.
    shortest_rr_interval_small_squares = 7
    duration_s = shortest_rr_interval_small_squares * 0.04
    rate_bpm = 60 / duration_s
    print(f"Observation: The rate is very fast. The shortest R-R interval is approximately {shortest_rr_interval_small_squares} small squares.")
    print(f"Calculation: Interval = {shortest_rr_interval_small_squares} * 0.04s = {duration_s:.2f}s")
    print(f"Peak Rate = 60 / {duration_s:.2f}s = {rate_bpm:.0f} bpm.")
    print("Conclusion: The heart rate is tachycardic, exceeding 200 bpm at times.\n")

    print("Step 3: QRS Complex Analysis")
    # A normal QRS is less than 3 small squares (0.12s). These are wider.
    qrs_width_small_squares = 4
    qrs_duration_s = qrs_width_small_squares * 0.04
    print(f"Observation: The QRS complexes are wide (e.g., >{qrs_width_small_squares-1} small squares or {qrs_duration_s-0.04:.2f}s) and their shape changes from beat to beat.")
    print("Conclusion: This indicates a wide-complex tachycardia with variable ventricular activation.\n")

    print("Step 4: Synthesis and Diagnosis")
    print("The combination of an irregularly irregular rhythm, a very fast rate (>200 bpm), and wide, variable QRS complexes is the classic presentation of Pre-excited Atrial Fibrillation.")
    print("This occurs in patients with an accessory pathway (like WPW syndrome) where atrial impulses are conducted rapidly to the ventricles, bypassing the AV node's rate-limiting effect.")
    print("\nFinal Diagnosis Choice: Pre-excited Atrial Fibrillation")

analyze_ecg()