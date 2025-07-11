def analyze_seizure_eeg():
    """
    Analyzes the EEG and EMG data to classify the seizure type.
    """
    # Step 1: Analyze EEG Findings
    eeg_morphology = "Generalized high-amplitude polyspike-and-wave discharges."
    eeg_frequency_hz = "Approximately 4-5 Hz."

    # Step 2: Analyze EMG Findings
    emg_pattern = "Brief, repetitive muscle jerks."
    eeg_emg_correlation = "EMG jerks are time-locked with the polyspikes in the EEG."

    # Step 3: Evaluate seizure types
    analysis_summary = """
The observed pattern is characteristic of a myoclonic seizure. Here's why:
- The EEG shows generalized polyspike-and-wave discharges, a common finding in myoclonic events.
- The EMG confirms the presence of brief, shock-like muscle contractions (myoclonus).
- The tight time-locking between the polyspikes (cortical discharge) and the EMG bursts (muscle jerks) is the hallmark of cortically-derived myoclonus.
- Other seizure types are ruled out:
  - No sustained muscle contraction for a 'tonic' seizure.
  - The motor component is too prominent for a typical 'absence' seizure.
  - There is no initial tonic phase, ruling out a 'tonic-clonic' seizure.
"""

    # Step 4: Final Answer
    final_answer = "C" # Myoclonic Seizure

    print("Seizure Analysis:")
    print("-----------------")
    print(f"EEG Morphology: {eeg_morphology}")
    print(f"EEG Frequency: {eeg_frequency_hz}")
    print(f"EMG Pattern: {emg_pattern}")
    print(f"EEG-EMG Correlation: {eeg_emg_correlation}")
    print("\nConclusion:")
    print(analysis_summary)
    print(f"Based on the evidence, the most accurate classification is 'Myoclonic Seizure'.")
    print("\nFinal Answer Code:")
    print(f'<<<{final_answer}>>>')

analyze_seizure_eeg()