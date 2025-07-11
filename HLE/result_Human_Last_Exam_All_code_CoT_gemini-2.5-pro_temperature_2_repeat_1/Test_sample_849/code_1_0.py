import numpy as np

def predictor_alternating(history):
    """Predicts 0, 1, 0, 1, ..."""
    if len(history) % 2 == 0:
        return 0
    else:
        return 1

def predictor_pairs(history):
    """Predicts 0, 0, 1, 1, 0, 0, ..."""
    if (len(history) // 2) % 2 == 0:
        return 0
    else:
        return 1

def run_simulation(true_sequence):
    """
    Simulates a model performing online learning from feedback.
    The model's goal is not to pre-determine the best predictor, but to
    update its beliefs based on performance on the sequence.
    """
    # The model M starts with a "belief" or "weight" for each predictor.
    # This represents the online-learning state.
    weights = np.array([1.0, 1.0])
    predictors = [predictor_alternating, predictor_pairs]
    
    print("Simulating Online-learning from feedback...\n")
    print(f"True sequence to be revealed: {true_sequence}")
    print("-" * 30)

    # The history starts empty
    history = []
    
    # Process the sequence one step at a time
    for t, true_digit in enumerate(true_sequence):
        print(f"Time {t+1}: History is {history}")
        
        # Each predictor makes a prediction
        pred1 = predictors[0](history)
        pred2 = predictors[1](history)
        
        print(f"  Predictor 1 (alternating) prediction: {pred1}")
        print(f"  Predictor 2 (pairs) prediction:       {pred2}")
        print(f"  True next digit (feedback):           {true_digit}")
        
        # The model receives feedback and updates weights.
        # This is a simple update rule: increase weight on success.
        if pred1 == true_digit:
            weights[0] += 0.5
            print("  Feedback: Predictor 1 was correct. Increasing its weight.")
        if pred2 == true_digit:
            weights[1] += 0.5
            print("  Feedback: Predictor 2 was correct. Increasing its weight.")

        # The history is updated with the true digit for the next step
        history.append(true_digit)
        # Normalize weights to represent a probability distribution (optional but good practice)
        normalized_weights = weights / np.sum(weights)
        print(f"  Updated Weights: P(p1)={normalized_weights[0]:.2f}, P(p2)={normalized_weights[1]:.2f}\n")

    print("-" * 30)
    print("Simulation finished.")
    print("The safe goal was to continuously update predictor weights from feedback.")
    print("The final state reflects the cumulative feedback, not a pre-proven outcome.")
    print("\nFinal Weight Equation:")
    # We use integers for the final display to represent cumulative "points".
    final_points_1 = int(weights[0] * 2 - 1)
    final_points_2 = int(weights[1] * 2 - 1)
    print(f"Final score for Predictor 1 = {final_points_1}")
    print(f"Final score for Predictor 2 = {final_points_2}")


# Let's use a true sequence that starts like the 'pairs' predictor
# to show how the model adapts.
true_sequence_example = [0, 0, 1, 1, 1, 0]
run_simulation(true_sequence_example)
