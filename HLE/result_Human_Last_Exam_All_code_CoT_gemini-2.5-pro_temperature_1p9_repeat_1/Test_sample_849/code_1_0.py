import collections

def simulate_conservative_learning():
    """
    Simulates a model performing conservative learning from a set of hypotheses.
    The model maintains a set of active predictors and discards those that are
    falsified by the observed data at each step.
    """
    # Define the "ground truth" sequence the model is trying to predict.
    # The initial part [1, 1, 0] is provided, and the rest must be predicted.
    target_sequence = [1, 1, 0, 1, 1, 0, 1, 0]

    # Define a set of hypotheses, represented by predictor functions.
    # Each function takes a time step `t` (index) and predicts the digit.
    predictors = {
        0: ("Always predicts 1", lambda t: 1),
        1: ("Alternates 1 and 0", lambda t: t % 2),
        2: ("Predicts '110' pattern", lambda t: [1, 1, 0][t % 3]),
        3: ("Always predicts 0", lambda t: 0),
        4: ("A complex but incorrect rule", lambda t: 1 if t in [0, 1, 3] else 0),
    }

    print("Starting simulation of Conservative Learning...")
    print(f"Target Sequence: {target_sequence}")
    print("Hypotheses (Predictors):")
    for idx, (desc, _) in predictors.items():
        print(f"  ID {idx}: {desc}")
    print("-" * 30)

    # The model starts with all hypotheses being active.
    # This set represents the numbers (indices) of the current valid hypotheses.
    active_hypotheses = set(predictors.keys())
    print(f"Initial set of active hypothesis indices: {sorted(list(active_hypotheses))}\n")

    for t, actual_digit in enumerate(target_sequence):
        if not active_hypotheses:
            print(f"Time {t}: No active hypotheses remain. Cannot predict.")
            break

        print(f"--- Time Step {t} ---")
        print(f"Observed Digit: {actual_digit}")
        print(f"Current active hypotheses: {sorted(list(active_hypotheses))}")

        falsified_hypotheses = set()

        # Evaluate each active hypothesis against the observed digit.
        for hypo_idx in active_hypotheses:
            _, predictor_func = predictors[hypo_idx]
            prediction = predictor_func(t)
            
            # This check represents the core "learning" step.
            if prediction != actual_digit:
                falsified_hypotheses.add(hypo_idx)

        # Update the set of active hypotheses by removing the falsified ones.
        if falsified_hypotheses:
            print(f"Falsified hypotheses at this step: {sorted(list(falsified_hypotheses))}")
            active_hypotheses -= falsified_hypotheses
        else:
            print("No hypotheses were falsified at this step.")
        
        # The final set for the next round. Each number is an index of a non-falsified hypothesis.
        print(f"Equation of state: Set of remaining valid hypotheses = {sorted(list(active_hypotheses))}\n")

    print("--- Simulation End ---")
    print(f"Final set of non-falsified hypotheses: {sorted(list(active_hypotheses))}")

simulate_conservative_learning()