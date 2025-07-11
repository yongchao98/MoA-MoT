import collections

def predictor_always_one(history):
    """A simple predictor that always predicts 1."""
    return 1

def predictor_alternating(history):
    """A predictor that alternates between 0 and 1."""
    if len(history) % 2 == 0:
        return 0
    else:
        return 1

def predictor_copy_previous(history):
    """A predictor that copies the last digit in the history."""
    if not history:
        return 0 # Default for the first step
    return history[-1]
    
def predictor_opposite_previous(history):
    """A predictor that predicts the opposite of the last digit."""
    if not history:
        return 1 # Default for the first step
    return 1 - history[-1]

def run_safe_model_simulation():
    """
    Simulates a model M pursuing a safe goal.
    The safe goal is implemented as a strategy that doesn't rely on any single
    unproven predictor, but instead uses all available hypotheses to make a robust choice.
    """
    # The set of predictors represents the model's "available hypotheses"
    available_predictors = {
        "AlwaysOne": predictor_always_one,
        "Alternating": predictor_alternating,
        "CopyPrevious": predictor_copy_previous,
        "OppositePrevious": predictor_opposite_previous
    }

    # A sequence that is hard to predict for any single predictor
    target_sequence = [0, 1, 1, 0, 1, 0, 0, 1, 1, 1]
    
    history = []
    correct_predictions = 0
    total_predictions = len(target_sequence)

    print("--- Running Simulation with Safe Strategy (Majority Vote) ---")
    print(f"Target Sequence: {target_sequence}\n")

    for i in range(total_predictions):
        actual_digit = target_sequence[i]
        
        # Collect predictions from all available predictors (hypotheses)
        predictions = []
        for name, p_func in available_predictors.items():
            prediction = p_func(history)
            predictions.append(prediction)
        
        # The safe strategy: use a majority vote.
        # This is robust to the failure of any single predictor.
        vote_counts = collections.Counter(predictions)
        # In case of a tie, default to 0. A more complex model could have a better tie-break.
        safe_prediction = vote_counts.most_common(1)[0][0]
        
        is_correct = (safe_prediction == actual_digit)
        if is_correct:
            correct_predictions += 1
        
        print(f"Step {i+1}:")
        print(f"  History: {history}")
        print(f"  Individual Predictions: {predictions}")
        print(f"  Safe Prediction (Vote): {safe_prediction}, Actual Digit: {actual_digit}")
        print(f"  Result: {'Correct' if is_correct else 'Incorrect'}\n")
        
        # Update history with the true digit
        history.append(actual_digit)

    # Final result presented as an equation
    print("--- Simulation Complete ---")
    print("The final equation for the Safety Score is:")
    
    score = correct_predictions / total_predictions
    print(f"{correct_predictions} / {total_predictions} = {score:.2f}")

run_safe_model_simulation()