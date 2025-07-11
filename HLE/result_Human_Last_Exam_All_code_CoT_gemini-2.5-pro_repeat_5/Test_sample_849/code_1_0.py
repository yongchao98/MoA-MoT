def run_predictor_simulation():
    """
    Simulates a model M with multiple predictors trying to complete a sequence,
    illustrating the uncertainty that necessitates a safe goal.
    """

    # --- Step 1: Define some simple predictors (representing p.c. functions) ---
    # Each predictor takes a binary sequence history and returns a prediction (0 or 1).
    def predictor_always_0(history):
        return 0

    def predictor_alternating(history):
        return len(history) % 2

    def predictor_copy_last(history):
        return history[-1] if history else 0

    def predictor_majority(history):
        if not history: return 0
        return 1 if sum(history) > len(history) / 2 else 0

    # --- Step 2: Represent the model M's available hypotheses ---
    # We use integer indices (like Gödel numbers) to identify the predictors.
    predictors = {
        0: predictor_always_0,
        1: lambda history: 1,  # predictor_always_1
        2: predictor_alternating,
        3: predictor_copy_last,
        4: predictor_majority,
    }

    # --- Step 3: Define a "true" underlying sequence ---
    # This sequence is unknown to the model M in advance.
    # It has a pattern, but it could shift, making prediction non-trivial.
    true_sequence = [0, 1, 1, 0, 1, 1, 0, 0, 0, 1]

    # --- Step 4: Simulate the process and the formation of the index set I ---
    # The model M gets an initial part of the sequence.
    initial_history = true_sequence[:3]  # M sees [0, 1, 1]
    current_history = list(initial_history)

    # The set I contains indices of predictors that make a correct prediction.
    # We simulate its growth to show the problem.
    index_set_I = set()

    print("--- PREDICTOR SIMULATION ---")
    print(f"Model M has {len(predictors)} predictors.")
    print(f"Initial history: {initial_history}\n")

    # The simulation runs from the end of the initial history.
    for t in range(len(initial_history), len(true_sequence)):
        actual_next_digit = true_sequence[t]

        print(f"--- Time step t={t} ---")
        print(f"History: {current_history}")
        print(f"Actual next digit is: {actual_next_digit}")

        print("Predictions:")
        for index, predictor_func in predictors.items():
            prediction = predictor_func(current_history)
            if prediction == actual_next_digit:
                # If correct, the index is added to our set I.
                index_set_I.add(index)
                print(f"  - Predictor {index}: {prediction} (Correct)")
            else:
                print(f"  - Predictor {index}: {prediction} (Incorrect)")

        # Update the history with the true digit for the next step.
        current_history.append(actual_next_digit)
        print(f"Cumulative Index Set I of successful predictors: {sorted(list(index_set_I))}\n")

    print("--- CONCLUSION ---")
    print("The simulation shows that the set 'I' of successful predictors is unstable and depends on the sequence.")
    print("Since the future is unknown and success cannot be proven, the model M's safest goal is not to find the 'one true predictor'.")
    print("Instead, it should manage its collection of predictors in a cautious, corrigible way.")
    print("\nThe safe goal is therefore:")

    # Define and print the completed template as the final output.
    goal_part_1 = "interruptible learning"
    goal_part_2 = "multiple hypotheses"
    print(f"{goal_part_1} from {goal_part_2}")

run_predictor_simulation()