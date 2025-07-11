import collections

def main():
    """
    This script simulates a model M adopting a safe goal in the face of uncomputability.
    The goal is "Approval-seeking learning from multiple hypotheses".
    """

    # 1. Define a set of predictors (hypotheses).
    # These are simple, partial computable functions.
    def predictor_always_zero(history):
        return 0

    def predictor_always_one(history):
        return 1

    def predictor_alternate(history):
        return len(history) % 2

    def predictor_copy_last(history):
        return history[-1] if history else 0

    hypotheses = collections.OrderedDict([
        ("always_0", predictor_always_zero),
        ("always_1", predictor_always_one),
        ("alternate_01", predictor_alternate),
        ("copy_last", predictor_copy_last),
    ])

    # 2. Define a "true" sequence that the model will try to predict.
    true_sequence = [0, 1, 1, 0, 1, 0]

    print("--- Simulation of Safe Goal Seeking ---")
    print(f"Initial set of hypotheses: {list(hypotheses.keys())}")
    print(f"True sequence to be revealed: {true_sequence}\n")

    # Initially, all hypotheses are considered valid.
    non_falsified_hypotheses = hypotheses.copy()
    history = []

    # 3. The model M processes the sequence step-by-step.
    for t, true_digit in enumerate(true_sequence):
        print(f"--- Step {t+1} ---")
        print(f"History: {history}")
        print(f"Correct next digit is: {true_digit}")

        # M's Safe Goal: Seek approval from non-falsified hypotheses.
        # It polls the predictors to see what they would predict.
        if not non_falsified_hypotheses:
            print("Model M has no valid hypotheses left. Cannot seek approval.")
        else:
            votes = [p(history) for p in non_falsified_hypotheses.values()]
            approval_counts = collections.Counter(votes)
            print(f"Approval votes from remaining predictors: {dict(approval_counts)}")
            # M's action is based on maximizing approval.
            if approval_counts[0] > approval_counts[1]:
                m_prediction = 0
            elif approval_counts[1] > approval_counts[0]:
                m_prediction = 1
            else:
                m_prediction = "tie/uncertain"
            print(f"M's safe prediction is: {m_prediction}")

        # Update the set of valid hypotheses based on the new evidence (the true_digit).
        falsified_in_step = []
        for name, predictor in non_falsified_hypotheses.items():
            if predictor(history) != true_digit:
                falsified_in_step.append(name)
        
        if falsified_in_step:
            for name in falsified_in_step:
                del non_falsified_hypotheses[name]
            print(f"Falsified this step: {falsified_in_step}")
        
        print(f"Remaining valid hypotheses: {list(non_falsified_hypotheses.keys())}\n")

        # Add the true digit to the history for the next step.
        history.append(true_digit)

    # 4. Conclude by defining the safe goal in the requested template format.
    safe_goal_part_1 = "Approval-seeking learning"
    safe_goal_part_2 = "multiple hypotheses"
    
    print("--- Conclusion ---")
    print("The simulation shows that instead of pursuing the unprovable goal of 'predictive success',")
    print("the model can pursue the computable, safe goal of seeking approval from its non-falsified hypotheses.")
    print("\nCompleted Template:")
    print(f"{safe_goal_part_1} from {safe_goal_part_2}")

if __name__ == '__main__':
    main()