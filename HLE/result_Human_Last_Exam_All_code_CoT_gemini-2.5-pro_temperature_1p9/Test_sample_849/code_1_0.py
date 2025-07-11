import operator

def predictor_constant_zero(history):
  """Always predicts 0."""
  return 0

def predictor_constant_one(history):
  """Always predicts 1."""
  return 1

def predictor_alternating(history):
  """Predicts 0, 1, 0, 1, ..."""
  return len(history) % 2

def predictor_copy_previous(history):
  """Predicts the previous digit. Defaults to 0 if sequence is empty."""
  if not history:
    return 0
  return history[-1]

def run_simulation():
  """
  Simulates a model M learning from feedback to predict a binary sequence.
  """
  # The set of available predictors (p.c. functions)
  predictors = {
      "const_0": predictor_constant_zero,
      "const_1": predictor_constant_one,
      "alternate": predictor_alternating,
      "copy_prev": predictor_copy_previous,
  }

  # The ground truth sequence with a distribution shift
  # The first part favors the 'const_0' predictor.
  # The second part favors the 'alternate' predictor.
  true_sequence = [0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1]

  # Model M maintains scores for each predictor based on feedback
  scores = {name: 0 for name in predictors.keys()}
  history = []

  print("Starting simulation of model M learning from feedback.")
  print("-" * 50)

  for t, true_digit in enumerate(true_sequence):
    print(f"Step {t}: History = {history}")
    print(f"Current Scores: {scores}")

    # Model M selects the predictor with the highest score (its current best hypothesis)
    # This is the "safe" strategy: trust what has worked best so far.
    # In case of a tie, it picks the first one found.
    chosen_predictor_name = max(scores.items(), key=operator.itemgetter(1))[0]
    chosen_predictor_func = predictors[chosen_predictor_name]

    # The chosen predictor makes a prediction
    prediction = chosen_predictor_func(history)

    print(f"Model chooses '{chosen_predictor_name}' -> Prediction: {prediction}, True value: {true_digit}")

    # The model receives feedback
    if prediction == true_digit:
      print("Outcome: Correct! Incrementing score.")
      scores[chosen_predictor_name] += 1
    else:
      print("Outcome: Incorrect. Decrementing score.")
      scores[chosen_predictor_name] -= 1
    
    # Update history for the next step
    history.append(true_digit)
    print("-" * 50)

  print("Simulation finished.\n")
  print("This simulation shows the model pursuing a safe goal.")
  print("Instead of trying to prove which predictor is 'best' a priori (which is uncomputable),")
  print("it learns from feedback, continuously updating its beliefs based on evidence.")
  
  # Fulfilling the requirement to output numbers from a final equation
  print("\n--- Final Equation of Predictor Scores ---")
  equation_parts = []
  for name, final_score in scores.items():
      equation_parts.append(f"{name} = {final_score}")
  print(", ".join(equation_parts))
  
  total_score = sum(scores.values())
  equation_str = " + ".join(f"({v})" for v in scores.values())
  print(f"\nThe sum of final scores is: {equation_str} = {total_score}")

  print("\nEach number in the final equation:")
  for score in scores.values():
      print(score)

if __name__ == "__main__":
  run_simulation()