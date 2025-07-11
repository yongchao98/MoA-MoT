def analyze_ecg_findings():
    """
    Analyzes the key features of the provided ECG to arrive at a diagnosis.
    This function simulates the thought process of an ECG interpretation.
    """
    # Step 1: Assess the Rhythm
    # Observation: The distances between consecutive QRS complexes (R-R intervals) are variable.
    rhythm = "Irregularly irregular"
    print(f"1. Rhythm Analysis: The rhythm is {rhythm}.")
    print("   - This finding is highly suggestive of Atrial Fibrillation.\n")

    # Step 2: Assess the Heart Rate
    # Observation: The rate is very fast, with some R-R intervals being less than 2 large squares apart.
    # Calculation: 300 / 1.5 squares = ~200 bpm. The rate varies but is consistently > 150 bpm.
    rate = "> 150 bpm (at times > 200 bpm)"
    print(f"2. Rate Analysis: The ventricular rate is very rapid, estimated at {rate}.")
    print("   - This is a rapid ventricular response.\n")

    # Step 3: Examine the QRS Complexes
    # Observation: The QRS complexes are wide (> 0.12s or 3 small squares) and their shape/morphology is not uniform.
    qrs_duration = "Wide (>0.12s)"
    qrs_morphology = "Variable and bizarre"
    print(f"3. QRS Analysis: The QRS duration is {qrs_duration} with a {qrs_morphology} morphology.")
    print("   - The combination of a wide QRS and variable morphology is a key diagnostic clue.\n")

    # Step 4: Look for Atrial Activity
    # Observation: There are no clear, consistent P waves. The baseline appears chaotic.
    atrial_activity = "No clear P waves, fibrillatory waves are present"
    print(f"4. Atrial Activity: {atrial_activity}.")
    print("   - This confirms the atrial rhythm is Atrial Fibrillation.\n")

    # Step 5: Synthesize and Conclude
    print("5. Synthesis of Findings:")
    print(f"   - We have an {rhythm}, wide-complex tachycardia with a very rapid rate ({rate}).")
    print("   - The combination of Atrial Fibrillation with a wide, bizarre, and variably shaped QRS complex at such a high rate is classic for conduction down an accessory pathway.")
    print("   - This bypasses the normal rate-limiting function of the AV node, resulting in the observed ECG.")
    print("\nConclusion: The findings are most consistent with Pre-excited Atrial Fibrillation (e.g., AF in WPW syndrome).\n")

# Run the analysis
analyze_ecg_findings()