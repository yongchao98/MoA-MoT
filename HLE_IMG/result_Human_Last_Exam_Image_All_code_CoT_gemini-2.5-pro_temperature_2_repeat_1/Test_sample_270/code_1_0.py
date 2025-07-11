import textwrap

def analyze_eeg_trace():
    """
    Analyzes the provided EEG and EMG traces to classify the seizure type.
    """
    analysis = """
    1. EEG Analysis (L-FC, L-SC, R-SC): The EEG traces show a sudden onset of high-amplitude, rhythmic discharges that are generalized (present in all three cortical channels simultaneously). The pattern consists of bursts of multiple spikes followed by a slow wave (polyspike-and-wave discharges). The frequency of these bursts is approximately 3-4 Hz (3-4 bursts per second, estimated from the 1s time bar).

    2. EMG Analysis: The EMG trace shows rhythmic bursts of electrical activity, indicating sharp, repetitive muscle contractions.

    3. Correlation of EEG and EMG: Crucially, each burst of muscle activity on the EMG is precisely time-locked with a polyspike-and-wave discharge on the EEG channels.

    4. Classification:
        - This pattern is not a Tonic seizure, as the EMG shows rhythmic bursts, not a sustained contraction.
        - It is not a typical Absence seizure, which usually lacks a significant motor component like the one seen here.
        - It is not a Tonic-clonic seizure, as there is no distinct initial tonic (stiffening) phase.
        - It perfectly matches the definition of a Myoclonic Seizure, which is characterized by brief, shock-like muscle jerks (myoclonus) that correspond directly to polyspike-and-wave discharges on the EEG. The brief, sharp nature of the EMG bursts reflects the myoclonic jerks.

    Therefore, the most accurate classification is Myoclonic Seizure.
    """
    
    print(textwrap.dedent(analysis).strip())
    print("\nFinal Answer:")

analyze_eeg_trace()

print("<<<C>>>")