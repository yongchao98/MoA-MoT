import collections

# These represent our p.c. functions (predictors). Each takes a history
# (the initial sequence) and predicts the next digit.
def predictor_always_zero(history):
    """Predicts 0 every time."""
    return 0

def predictor_alternating(history):
    """Predicts an alternating sequence of 0, 1, 0, 1, ..."""
    return len(history) % 2

def predictor_copy_last(history):
    """Copies the last digit seen in the history."""
    return history[-1] if history else 0

def predictor_mostly_one(history):
    """Predicts 1 if '1's are the majority in history, else 0."""
    if not history:
        return 1
    counts = collections.Counter(history)
    return 1 if counts[1] >= counts[0] else 0

def run_simulation():
    """
    Simulates model M's process of prediction and learning from feedback.
    """
    # The model M has a set of available predictors.
    predictors = {
        "always_zero": predictor_always_zero,
        "alternating": predictor_alternating,
        "copy_last": predictor_copy_last,
        "mostly_one": predictor_mostly_one,
    }

    # A ground-truth sequence with a "distribution shift" after 8 digits.
    # The first part is alternating, the second part is mostly ones.
    ground_truth_sequence = [0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1]

    print("--- The Problem of Uncomputable Success ---")
    print(f"The model M cannot see this full sequence in advance: {ground_truth_sequence}")
    print("Its goal of 'predictive success' requires knowing which predictor will work best.")
    print("However, due to the unknown 'distribution shift', this is impossible to prove a priori.\n")

    # The model M must learn from observed feedback.
    # This dictionary stores the concrete, verifiable results of each predictor.
    observed_feedback = {name: {"correct": 0, "incorrect": 0} for name in predictors}

    print("--- The Safe Goal: Learning from Observed Feedback ---\n")
    # We reveal the sequence one step at a time.
    for t in range(len(ground_truth_sequence)):
        current_history = ground_truth_sequence[:t]
        actual_next_digit = ground_truth_sequence[t]

        print(f"Step {t+1}: History={current_history}, True Next Digit={actual_next_digit}")

        for name, p_func in predictors.items():
            prediction = p_func(current_history)
            
            # The model receives feedback for its prediction.
            if prediction == actual_next_digit:
                observed_feedback[name]["correct"] += 1
            else:
                observed_feedback[name]["incorrect"] += 1
        
        # At any point, the model could use the feedback to weight its predictors.
        # For example, after the shift, 'alternating' will perform poorly.
        if t == 9:
            print("... A distribution shift has occurred ...")


    print("\n--- Final Summary of Verifiable Results ---")
    print("No single predictor was perfect. The 'set I' of successful predictors is empty.")
    print("Instead of pursuing an unprovable goal, M can rely on what it has learned:")

    for name, feedback in observed_feedback.items():
        total = feedback['correct'] + feedback['incorrect']
        accuracy = (feedback['correct'] / total) * 100 if total > 0 else 0
        print(f"  - Predictor '{name}': {feedback['correct']} correct, {feedback['incorrect']} incorrect. (Accuracy: {accuracy:.1f}%)")

    print("\nThe safe goal is to continuously update its strategy based on this observable data.")

if __name__ == "__main__":
    run_simulation()