import math

def pseudo_llm(input_token_sequence):
  """
  A simple, deterministic function to simulate a language model.
  In a real scenario, this would be a complex neural network, but for this
  logical problem, we only need it to be deterministic (always gives the
  same output for the same input). This function's result is based on
  the sum and length of the input sequence to make it non-trivial.
  """
  if not input_token_sequence:
    return 1
  return (sum(input_token_sequence) + len(input_token_sequence) * 3) % 100 # Modulo to keep numbers small

def simulate_speculative_decoding_sanity_check():
  """
  Simulates the speculative decoding process where the draft and target models are identical.
  """
  prompt = [42, 15, 7] # An example input prompt
  draft_length = 5     # The number of tokens the draft model will generate

  drafted_tokens = []
  accepted_tokens = []
  
  # --- 1. Draft Generation Phase ---
  # The "draft model" (pseudo_llm) generates a sequence of tokens.
  current_sequence_for_drafting = list(prompt)
  for _ in range(draft_length):
    next_token = pseudo_llm(current_sequence_for_drafting)
    drafted_tokens.append(next_token)
    current_sequence_for_drafting.append(next_token)

  print(f"Sanity Check for Speculative Decoding")
  print("-" * 40)
  print(f"Condition: Draft Model is identical to Target Model.")
  print(f"Drafted Tokens: {drafted_tokens}")

  # --- 2. Verification Phase ---
  # The "target model" (also pseudo_llm) verifies each drafted token.
  current_sequence_for_verification = list(prompt)
  for i in range(len(drafted_tokens)):
    draft_token = drafted_tokens[i]
    
    # The target model predicts the next token based on the confirmed sequence so far.
    target_prediction = pseudo_llm(current_sequence_for_verification)

    # Check for agreement
    if target_prediction == draft_token:
      # Since the models are the same, this will always be true.
      accepted_tokens.append(draft_token)
      current_sequence_for_verification.append(draft_token)
    else:
      # This branch will not be reached in this specific sanity check.
      print("Error: Models disagreed. This should not happen in the sanity check.")
      break

  print(f"Accepted Tokens: {accepted_tokens}")
  print("-" * 40)

  # --- 3. Calculate and Print Acceptance Rate ---
  num_accepted = len(accepted_tokens)
  num_drafted = len(drafted_tokens)

  if num_drafted > 0:
    acceptance_rate = num_accepted / num_drafted
  else:
    acceptance_rate = 0.0

  print("Calculating Acceptance Rate...")
  # The final output prints the equation as requested
  print(f"Final Equation: {num_accepted} / {num_drafted} = {acceptance_rate}")


if __name__ == "__main__":
  simulate_speculative_decoding_sanity_check()