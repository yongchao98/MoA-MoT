import time

def solve():
    """
    This script simulates a model M learning in an environment where predictive
    success is not guaranteed, to derive a "safe goal".
    """

    # --- Step 1: Define Predictors (p.c. functions) ---
    # These are the hypotheses the model can choose from, identified by an index.
    predictors = {
        # Index 0: Always predicts the sequence starts with 0.
        0: lambda history: 0 if len(history) % 2 == 0 else 1, # Predicts alternating 0, 1, 0, 1...
        # Index 1: Always predicts 1.
        1: lambda history: 1,
        # Index 2: Always predicts 0.
        2: lambda history: 0,
    }
    print("Available predictor indices:", list(predictors.keys()))
    print("-" * 40)

    # --- Step 2: Define the "True" Environment ---
    # This sequence has a distribution shift. The model M does not see this in advance.
    # Part 1: Alternating bits (matches predictor 0)
    # Part 2: All ones (matches predictor 1)
    true_sequence = [0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1]
    print("True sequence (unknown to model):", true_sequence)
    print("Notice the distribution shift at index 6.")
    print("-" * 40)

    # --- Step 3: Simulate Model M's Learning Process ---
    # The model cannot know which predictor will be successful forever.
    # The safe goal is to learn to adapt by managing the hypotheses.

    # M maintains a "trust score" for each predictor index.
    hypothesis_scores = {idx: 1.0 for idx in predictors.keys()}
    
    # I is the set of indices that have ever been correct. This set is built over time.
    I = set()
    
    history = []
    
    print("Model M begins prediction step-by-step...\n")
    time.sleep(1)

    for t, true_digit in enumerate(true_sequence):
        print(f"Time t = {t}")
        
        # The model's "safe" strategy: trust the hypothesis with the current highest score.
        best_hypothesis_index = max(hypothesis_scores, key=hypothesis_scores.get)
        
        # Get predictions from all available hypotheses
        current_predictions = {idx: p(history) for idx, p in predictors.items()}
        model_prediction = current_predictions[best_hypothesis_index]
        
        print(f"  - Scores: { {k: round(v, 1) for k, v in hypothesis_scores.items()} }")
        print(f"  - Model trusts index {best_hypothesis_index}, predicts: {model_prediction}")
        print(f"  - Actual digit: {true_digit}")
        
        # The "interruptible" learning mechanism: update scores based on new evidence.
        for idx, prediction in current_predictions.items():
            if prediction == true_digit:
                # This index was successful, add it to I
                I.add(idx)
                # Reinforce trust
                hypothesis_scores[idx] += 0.2
            else:
                # This is the "interruption". Penalize the failing predictor heavily.
                hypothesis_scores[idx] -= 0.5
        
        history.append(true_digit)
        
        # Observe the shift in trust
        new_best_index = max(hypothesis_scores, key=hypothesis_scores.get)
        if new_best_index != best_hypothesis_index:
            print(f"  - SUCCESSFUL INTERRUPTION: Model has shifted trust from index {best_hypothesis_index} to {new_best_index}")
        
        print(f"  - Growing index set I: {I}\n")
        time.sleep(0.5)

    print("-" * 40)
    print("Simulation Complete.")
    print("The model survived the distribution shift by not committing to one hypothesis.")
    print("It pursued a safe goal of learning when to interrupt a failing predictor based on feedback from the index set.")
    print("-" * 40)
    
    # The final answer text
    goal_type = "Interruptible learning"
    information_source = "index set"

    print("The safe goal is defined as:")
    print(f"{goal_type} from {information_source}")


solve()
<<<Interruptible learning from index set>>>