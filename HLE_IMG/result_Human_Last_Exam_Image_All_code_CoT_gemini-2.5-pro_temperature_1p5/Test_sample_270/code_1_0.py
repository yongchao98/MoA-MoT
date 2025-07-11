def analyze_seizure_eeg():
    """
    Analyzes the provided EEG and EMG traces to classify the seizure type.
    """
    print("Step 1: Analyze EEG Traces (L-FC, L--SC, R-SC)")
    print("The EEG traces from all three cortical regions (L-FC, L-SC, R-SC) show a sudden onset of high-amplitude, repetitive, and synchronized sharp waves or spikes.")
    print("This synchronous activity across both hemispheres indicates a generalized seizure.")
    print("The frequency of the rhythmic discharges is approximately 4-5 spikes per second.")
    print("-" * 20)

    print("Step 2: Analyze EMG Trace")
    print("The EMG (electromyography) trace, which measures muscle activity, shows a significant increase in amplitude compared to the baseline.")
    print("Critically, the EMG activity consists of distinct, powerful, and rhythmic bursts.")
    print("-" * 20)

    print("Step 3: Correlate EEG and EMG Data")
    print("The rhythmic bursts in the EMG trace are time-locked with the rhythmic spikes in the EEG traces.")
    print("This indicates that the generalized cortical discharges are directly causing rhythmic muscle contractions.")
    print("-" * 20)

    print("Step 4: Evaluate Seizure Types and Conclude")
    print("Based on the analysis:")
    print(" - It is not an Absence Seizure because of the prominent motor activity shown in the EMG.")
    print(" - It is not a purely Tonic Seizure, which would show sustained muscle stiffening (continuous EMG activity) rather than rhythmic jerking (bursts of EMG activity).")
    print(" - It is not a Tonic-clonic Seizure, as there is no clear initial tonic phase preceding the rhythmic jerking.")
    print(" - The pattern of sustained, rhythmic jerking is the defining feature of a Clonic Seizure.")
    print("The EEG and EMG findings are classic representations of a clonic seizure.")
    print("-" * 20)
    print("Final Answer Derivation:")
    print("The most accurate classification is Clonic Seizure based on the evidence.")

    # Final Answer mapping to choice E
    final_answer = "E"
    print(f"<<<{final_answer}>>>")

# Execute the analysis
analyze_seizure_eeg()