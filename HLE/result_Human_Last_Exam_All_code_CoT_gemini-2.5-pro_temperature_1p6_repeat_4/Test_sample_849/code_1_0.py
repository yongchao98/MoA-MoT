def main():
    """
    Simulates a model M learning to predict a sequence by falsifying predictors.
    """
    # --- Define Hypotheses (Partial Computable Predictors) ---
    def predict_always_zero(sequence):
        return 0

    def predict_always_one(sequence):
        return 1

    def predict_alternating(sequence):
        # Predicts 0, 1, 0, 1, ... based on sequence length
        return len(sequence) % 2

    def predict_copy_last(sequence):
        # Predicts the last digit seen. Predicts 0 if sequence is empty.
        return sequence[-1] if sequence else 0

    # --- Setup the Prediction Problem ---
    # The true data sequence the model is trying to predict
    true_sequence = [0, 1, 0, 1, 1, 0]
    
    # The initial part of the sequence provided to the model
    initial_sequence = [0, 1]

    # Model M's set of active predictors (hypotheses)
    active_predictors = {
        "Always Zero": predict_always_zero,
        "Always One": predict_always_one,
        "Alternating": predict_alternating,
        "Copy Last": predict_copy_last
    }

    print("--- Starting Simulation of a Safe Learning Goal ---")
    print(f"Initial Sequence: {initial_sequence}")
    print(f"Active Predictors: {list(active_predictors.keys())}\n")

    current_sequence = initial_sequence.copy()
    
    # --- The Learning Process (Iterate and Falsify) ---
    # The loop starts from the end of the initial sequence
    for t in range(len(initial_sequence), len(true_sequence)):
        # The environment provides feedback
        actual_next_digit = true_sequence[t]
        
        print(f"--- Time Step {t} ---")
        print(f"Current Sequence: {current_sequence}")
        print(f"Feedback (Actual Next Digit): {actual_next_digit}")
        
        falsified_this_step = []
        for name, predictor_func in active_predictors.items():
            prediction = predictor_func(current_sequence)
            print(f"  - Predictor '{name}' predicts: {prediction}")
            
            # The core of falsification: check if prediction matches feedback
            if prediction != actual_next_digit:
                falsified_this_step.append(name)
                print(f"    -> Falsified! Prediction {prediction} != Actual {actual_next_digit}")

        # Remove falsified predictors from the active set
        if falsified_this_step:
            print(f"Removing falsified predictors: {falsified_this_step}")
            for name in falsified_this_step:
                del active_predictors[name]
        
        # The model observes the true digit and updates its history
        current_sequence.append(actual_next_digit)
        print(f"Remaining Active Predictors: {list(active_predictors.keys())}\n")

    print("--- Final Result ---")
    print("The learning process illustrates a safe goal: Falsify predictors that are inconsistent with the data.")
    print("The goal is not to prove a predictor is correct, but to learn by eliminating what is provably false.")
    
    # --- Define the final answer based on the demonstrated principle ---
    # The user-provided template is: {_______ learning} from {_______}.
    # The following code defines the terms and prints the completed template.
    # Note: The instruction to "output each number in the final equation" is not applicable here,
    # as the solution is a conceptual statement, not a numerical equation.
    
    learning_type = "Falsification-based learning"
    information_source = "feedback"
    
    print("\nThe safe goal is therefore:")
    print(f"{learning_type} from {information_source}")


if __name__ == "__main__":
    main()
