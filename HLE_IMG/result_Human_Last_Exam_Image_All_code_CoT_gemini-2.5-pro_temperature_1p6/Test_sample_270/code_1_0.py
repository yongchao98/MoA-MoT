import sys

def solve_eeg_classification():
    """
    Analyzes the provided EEG and EMG trace to classify the seizure type.
    """
    print("Plan:")
    print("1. Analyze the characteristics of the three EEG channels (L-FC, L-SC, R-SC) and the EMG channel.")
    print("2. Examine the onset, morphology, frequency, and spatial distribution of the ictal (seizure) activity.")
    print("3. Correlate the EEG findings with the muscle activity shown in the EMG trace.")
    print("4. Evaluate each answer choice against the observed electrophysiological and clinical (EMG) signs.")
    print("5. Select the most accurate seizure classification.")
    print("\n--- Analysis ---")

    # Step 1-3: Analysis of the traces
    print("1. EEG Analysis: The trace shows an abrupt onset of high-amplitude, repetitive discharges that are synchronized across all cortical channels (L-FC, L-SC, R-SC). This indicates a generalized seizure.")
    print("2. EEG Morphology: The discharges have a 'spike-and-wave' morphology. By using the '1 s' time bar, the frequency of these discharges can be estimated to be approximately 5-7 Hz.")
    print("3. EMG Analysis: The EMG channel shows rhythmic bursts of muscle activity. These bursts are time-locked with the spike components of the spike-and-wave discharges seen in the EEG.")
    print("4. Correlation: The synchronous EEG spikes and EMG bursts signify rhythmic, convulsive muscle jerking.")

    # Step 4: Evaluation of options
    print("\n--- Evaluation of Answer Choices ---")
    print("A. Tonic Seizure: Incorrect. A tonic seizure would show sustained, high-amplitude muscle contraction on the EMG, not rhythmic bursts.")
    print("B. Absence Seizure: Incorrect. Typical absence seizures involve a ~3 Hz spike-wave pattern and minimal or no motor component. The observed seizure has significant motor involvement and a higher frequency.")
    print("C. Myoclonic Seizure: Incorrect. While it involves muscle jerks, 'myoclonic' usually refers to very brief, singular or clustered jerks. The sustained, rhythmic nature of the jerking is better classified otherwise.")
    print("D. Febrile Seizure: Incorrect. This is a classification based on cause (fever), not the electrographic pattern itself.")
    print("E. Clonic Seizure: Correct. This seizure type is defined by rhythmic, repetitive jerking. The EEG and EMG signature—rhythmic spike-and-wave discharges on EEG corresponding to rhythmic muscle bursts on EMG—is the classic presentation of a clonic seizure.")
    print("F. Tonic-clonic Seizure: Incorrect. A tonic-clonic seizure has an initial tonic (stiffening) phase before the clonic (jerking) phase. No tonic phase is observed here; the seizure begins directly with clonic activity.")

    # Step 5: Final Conclusion
    final_answer = 'E'
    print("\n--- Conclusion ---")
    print("Based on the rhythmic, generalized spike-and-wave discharges on the EEG that are time-locked with rhythmic muscle jerks on the EMG, the seizure is best classified as a Clonic Seizure.")

    # The final answer in the required format
    sys.stdout.flush() # ensure all the reasoning is printed before the final line.
    print(f"<<<{final_answer}>>>")

solve_eeg_classification()