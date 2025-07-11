import numpy as np

def simulate_online_learning():
    """
    Simulates a model M learning to trust predictors based on interactive feedback.
    """

    # 1. Define the ground truth sequence and the initial part M sees.
    # The pattern is mostly '1010...', but with some shifts/noise.
    ground_truth_sequence = "101011010100101011101"
    initial_part = "101"
    completion_sequence = ground_truth_sequence[len(initial_part):]

    # 2. Define a set of predictors (p.c. functions).
    # Their internal logic is unknown to M. M only sees their predictions.
    predictors = {
        "alternates_10": lambda seq: '0' if seq[-1] == '1' else '1',
        "repeats_last": lambda seq: seq[-1],
        "always_one": lambda seq: '1',
        "opposite_of_start": lambda seq: '0' if seq[0] == '1' else '1'
    }

    # M's confidence scores for each predictor. It starts with equal trust.
    # This represents M's internal state. The goal is to learn from feedback
    # to make these scores reflect the true utility of the predictors.
    confidence_scores = {name: 1.0 for name in predictors}

    print("Model M starts with an initial sequence part: '{}'".format(initial_part))
    print("M must predict the completion: '{}'".format(completion_sequence))
    print("-" * 30)
    print("Initial confidence scores: {}".format(confidence_scores))
    print("-" * 30)


    current_sequence = list(initial_part)
    # 3. M performs online learning from interactive feedback.
    for i, true_digit in enumerate(completion_sequence):
        print(f"Step {i+1}: Current sequence is {''.join(current_sequence)}")
        print(f"M needs to predict the next digit. The true next digit is: {true_digit}")

        # Each predictor makes a prediction based on the current sequence.
        predictions = {name: func("".join(current_sequence)) for name, func in predictors.items()}
        print(f"  > Predictor outputs: {predictions}")

        # M gets interactive feedback by comparing predictions to the true digit.
        # It then updates its confidence scores online.
        for name, pred_digit in predictions.items():
            if pred_digit == true_digit:
                # Reward correct predictor
                confidence_scores[name] *= 1.2
                print(f"    - '{name}' was CORRECT. Confidence increased.")
            else:
                # Penalize incorrect predictor
                confidence_scores[name] *= 0.8
                print(f"    - '{name}' was WRONG. Confidence decreased.")
        
        # Add the true digit to the sequence for the next step.
        current_sequence.append(true_digit)
        
        # Normalize scores for clarity
        total_score = sum(confidence_scores.values())
        normalized_scores = {name: score / total_score for name, score in confidence_scores.items()}

        print(f"  > Updated normalized confidence: " + \
              ", ".join([f"{name}: {score:.2f}" for name, score in normalized_scores.items()]))
        print("-" * 30)
    
    print("Final learned confidence scores after processing the entire sequence:")
    final_total_score = sum(confidence_scores.values())
    final_normalized_scores = {name: score / final_total_score for name, score in confidence_scores.items()}
    
    # This final state is the result of the safe learning goal.
    # The model hasn't proven anything, but it has learned which predictor to trust.
    # We print the "final equation" which we interpret as the final state of the learned scores.
    for name, score in final_normalized_scores.items():
        print(f"Final Belief('{name}') = {score:.4f}")

# Execute the simulation
simulate_online_learning()