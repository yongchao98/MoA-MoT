import math

def analyze_ecg():
    """
    Analyzes the provided ECG and explains the diagnosis.
    """
    print("ECG Analysis Steps:")
    print("===================")

    # 1. Rhythm Analysis
    rhythm = "Irregularly irregular"
    print(f"1. Rhythm: The R-R intervals are highly variable. The rhythm is {rhythm}.")
    print("   This strongly suggests an underlying atrial fibrillation.")

    # 2. Rate Calculation
    # Shortest R-R interval is approx 1.2 large squares (1.2 * 0.2s = 0.24s)
    # Longest R-R interval is approx 2 large squares (2 * 0.2s = 0.4s)
    shortest_rr_interval_sec = 1.2 * 0.2
    longest_rr_interval_sec = 2.0 * 0.2
    max_rate = math.ceil(60 / shortest_rr_interval_sec)
    min_rate = math.floor(60 / longest_rr_interval_sec)
    print("\n2. Heart Rate: The rate is very fast and variable.")
    print(f"   - Shortest R-R interval is ~1.2 large squares.")
    print(f"   - Corresponding Maximum Rate: 60 / (1.2 * 0.2) = {max_rate} bpm.")
    print(f"   - Longest R-R interval is ~2.0 large squares.")
    print(f"   - Corresponding Minimum Rate: 60 / (2.0 * 0.2) = {min_rate} bpm.")
    print(f"   The estimated rate ranges from {min_rate} to {max_rate} bpm.")


    # 3. QRS Analysis
    qrs_duration_ms = 140 # Estimated average
    print("\n3. QRS Complex:")
    print(f"   - The QRS complexes are wide (>{qrs_duration_ms-20}ms).")
    print("   - The morphology (shape) of the QRS complexes is bizarre and varies from beat to beat.")

    # 4. Synthesize and Conclude
    print("\n4. Synthesis of Findings:")
    print("- Irregularly irregular rhythm -> Atrial Fibrillation.")
    print("- Wide QRS complexes -> Ventricular origin, aberrancy, or pre-excitation.")
    print("- Extremely fast rate (>200 bpm) AND varying QRS morphology.")
    print("\nThis specific combination of an irregularly irregular, extremely fast, wide-complex tachycardia with beat-to-beat morphological changes is the hallmark of Atrial Fibrillation with an accessory pathway (pre-excitation). The accessory pathway bypasses the AV node, allowing for dangerously rapid ventricular response.")
    
    print("\nConclusion:")
    print("The diagnosis is Pre-excited Atrial Fibrillation.")

if __name__ == "__main__":
    analyze_ecg()