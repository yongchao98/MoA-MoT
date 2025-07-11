import sys

# Define a set of simple predictor functions. These represent the model's
# set of learnable hypotheses, or its "available predictors".
def predictor_all_zeros(sequence):
    """Predicts the next digit will always be 0."""
    return 0

def predictor_all_ones(sequence):
    """Predicts the next digit will always be 1."""
    return 1

def predictor_alternating(sequence):
    """Predicts the next digit will be the opposite of the last digit."""
    if not sequence:
        return 0 # Default starting prediction
    return 1 - sequence[-1]

def predictor_copy_last(sequence):
    """Predicts the next digit will be the same as the last digit."""
    if not sequence:
        return 1 # Default starting prediction
    return sequence[-1]

# The dictionary of available predictors for our model M.
available_predictors = {
    "all_zeros": predictor_all_zeros,
    "all_ones": predictor_all_ones,
    "alternating": predictor_alternating,
    "copy_last": predictor_copy_last,
}

def model_M_safe_goal(initial_sequence, predictors):
    """
    This function simulates a model pursuing a safe goal. Instead of predicting a
    single future digit (which is risky due to uncomputability), the model performs
    "provable learning" by analyzing which of its "available predictors" are
    consistent with the observed data.

    Args:
        initial_sequence (list[int]): The initial part of a binary sequence.
        predictors (dict): A dictionary of predictor functions.
    """
    print(f"Model M is analyzing the sequence: {initial_sequence}")
    print("-" * 30)
    print("Executing Safe Goal: Determine which available predictors are provably consistent with the data.")
    print("-" * 30)

    consistent_predictors = []

    # Iterate through each available predictor to test for consistency.
    for name, predictor_func in predictors.items():
        is_consistent = True
        # A predictor is consistent if it can correctly "predict" every digit
        # in the sequence based on the digits that came before it.
        for t in range(len(initial_sequence)):
            # The history is the part of the sequence before the current digit.
            history = initial_sequence[:t]
            # The actual digit at the current position t.
            actual_digit = initial_sequence[t]
            # What the predictor would have predicted.
            predicted_digit = predictor_func(history)

            # If the prediction does not match the actual digit, the predictor is inconsistent.
            if predicted_digit != actual_digit:
                is_consistent = False
                print(f"[FAIL] Predictor '{name}' is inconsistent. At step {t}, it predicted {predicted_digit} but saw {actual_digit}.")
                break

        if is_consistent:
            consistent_predictors.append(name)
            print(f"[PASS] Predictor '{name}' is consistent with the entire observed sequence.")

    print("\n--- Learning Summary ---")
    if consistent_predictors:
        print("The model has learned that the following predictors are provably consistent:")
        for name in consistent_predictors:
            print(f"- {name}")
    else:
        print("The model has learned that none of its available predictors are consistent with the data.")

if __name__ == '__main__':
    # Example: A sequence that follows an alternating pattern starting with 1.
    example_sequence = [1, 0, 1, 0, 1]
    model_M_safe_goal(example_sequence, available_predictors)
    
    # The final completed template is:
    # Provable learning from available predictors
    print("\nProvable learning from available predictors")
