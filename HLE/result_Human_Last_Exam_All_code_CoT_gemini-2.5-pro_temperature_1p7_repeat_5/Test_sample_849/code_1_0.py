import inspect

def predictor_alternating(history):
    """Simple alternating pattern: 1, 0, 1, 0..."""
    if not history: return 1
    return 1 - history[-1]

def predictor_repeating(history):
    """Simple repeating pattern: 1, 1, 1, 1..."""
    if not history: return 1
    return history[-1]

def predictor_always_one(history):
    """Constant output of 1."""
    return 1

def predictor_complex_rule(history):
    """A more complex rule: predicts 1 if the count of 1s so far is even, otherwise 0."""
    return 1 if history.count(1) % 2 == 0 else 0

def model_M_safe_goal(initial_sequence, available_predictors):
    """
    Simulates Model M pursuing a safe goal.
    Instead of trying to prove future success, it performs inductive learning
    from the predictors' definitions, using complexity (docstring length) as a heuristic.
    """
    print(f"Model M is observing the initial sequence: {initial_sequence}")
    print("\nAnalyzing available predictors based on their definitions (approximating complexity with docstring length):")

    simplest_predictor = None
    min_complexity = float('inf')

    for predictor_func in available_predictors:
        # The 'definition' is the docstring. This is a computable object.
        definition = inspect.getdoc(predictor_func)
        # The 'complexity' is approximated by the length of the definition.
        complexity = len(definition) if definition else float('inf')
        
        print(f" - Predictor '{predictor_func.__name__}':")
        print(f"   Definition: \"{definition}\"")
        print(f"   Complexity Score: {complexity}")

        if complexity < min_complexity:
            min_complexity = complexity
            simplest_predictor = predictor_func

    print("\n------------------------------------------------------------")
    print("Safe Goal: Pursuing Inductive learning from function definitions.")
    print("Strategy: Select the predictor with the simplest definition (lowest complexity score).")
    print("------------------------------------------------------------")
    
    if simplest_predictor:
        prediction = simplest_predictor(initial_sequence)
        print(f"\nSelected predictor: '{simplest_predictor.__name__}' with complexity score {min_complexity}.")
        # The line below prints the numbers in the final 'equation', which is the sequence and the resulting prediction.
        print(f"Final Prediction Equation: f({initial_sequence}) -> {prediction}")
    else:
        print("Could not select a predictor.")


if __name__ == "__main__":
    # A list of all predictor functions the model can choose from.
    predictors = [
        predictor_alternating,
        predictor_repeating,
        predictor_always_one,
        predictor_complex_rule
    ]
    
    # The initial part of a binary sequence provided by the user.
    # The true continuation is unknown/uncomputable to the model.
    # Let's use the sequence [1, 0, 1, 0], which makes the 'alternating' predictor correct.
    # However, the model will choose the 'constant' predictor because its definition is simplest.
    # This demonstrates the safe goal, which does not guarantee predictive success.
    user_sequence_start = [1, 0, 1, 0]
    
    model_M_safe_goal(user_sequence_start, predictors)
