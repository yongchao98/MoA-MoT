def predictor_alternating(history):
    """Predicts alternating bits: 0, 1, 0, 1, ..."""
    if not history or len(history) % 2 == 0:
        return 0
    else:
        return 1

def predictor_always_one(history):
    """Always predicts 1."""
    return 1

def predictor_repeater(history):
    """Predicts the last bit, or 0 if history is empty."""
    if not history:
        return 0
    return history[-1]
    
def predictor_opposite(history):
    """Predicts the opposite of the last bit, or 1 if history is empty."""
    if not history:
        return 1
    return 1 - history[-1]

def run_simulation():
    """
    Simulates a model M pursuing a safe goal by filtering predictors based on evidence.
    """
    # The model M is given a set of hypotheses (predictors).
    # These are our "formal systems", identified by name.
    hypotheses = {
        "alternating": predictor_alternating,
        "always_one": predictor_always_one,
        "repeater": predictor_repeater,
        "opposite": predictor_opposite,
    }

    # The true sequence, unknown to the model, exhibits a distribution shift.
    # It starts by alternating (0,1,0,1) then switches to repeating the last bit (1,1,1...).
    true_sequence = [0, 1, 0, 1, 1, 1, 1] 

    # The model maintains a set of hypotheses consistent with evidence so far.
    # This is the practical aspect of the "provably-safe" goal.
    consistent_hypotheses = set(hypotheses.keys())

    history = []

    print("--- Simulation of Provably-Safe Learning from Formal Systems ---")
    print(f"Initial set of hypotheses: {sorted(list(consistent_hypotheses))}\n")

    # The model observes the sequence and updates its set of consistent hypotheses.
    for t, true_bit in enumerate(true_sequence):
        print(f"--- Time step t={t} ---")
        print(f"History: {history}")
        print(f"Observing next bit: {true_bit}")
        
        if not consistent_hypotheses:
            print("No remaining consistent hypotheses. The model's initial set was insufficient.")
            break
            
        falsified_this_step = set()
        
        # Check each currently-believed hypothesis against the new evidence.
        for name in consistent_hypotheses:
            predictor_func = hypotheses[name]
            prediction = predictor_func(history)
            
            # This check is a provable step: the prediction is either right or wrong.
            if prediction != true_bit:
                falsified_this_step.add(name)
        
        # The core of the safe learning process: remove what has been proven wrong.
        if falsified_this_step:
            print(f"Falsified hypotheses: {sorted(list(falsified_this_step))}")
            consistent_hypotheses.difference_update(falsified_this_step)
        
        print(f"Updated set of consistent hypotheses: {sorted(list(consistent_hypotheses))}")
        print("")

        history.append(true_bit)

    print("--- Final Result ---")
    print(f"The final sequence was: {history}")
    print(f"The final set of non-falsified hypotheses is: {sorted(list(consistent_hypotheses))}")

run_simulation()