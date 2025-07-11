import numpy as np

# In speculative decoding, a draft token is accepted based on the ratio of
# probabilities from the target and draft models. Let's simulate this.

# We need a mock model function that returns probabilities for next tokens.
# In a real scenario, this would be a complex neural network. For this
# demonstration, a simple dictionary lookup is sufficient to prove the concept.
def mock_language_model(context):
    """
    A mock Language Model. For a given input context, it returns a
    fixed probability distribution over a small vocabulary.
    This simulates a deterministic model's output for a given input.
    """
    # This distribution is arbitrary and just for demonstration.
    return {"the": 0.5, "a": 0.2, "cat": 0.1, "dog": 0.1, "runs": 0.1}

# The problem states that for our sanity check, the draft model and the
# target model are the same.
draft_model = mock_language_model
target_model = mock_language_model

# Let's assume the draft model generated a sequence of 3 tokens.
# The actual tokens don't matter as long as they are in the vocabulary.
drafted_tokens = ["the", "dog", "runs"]
context = "start_of_sentence"  # Context for the generation step

# We will count how many tokens are drafted and how many are accepted.
num_drafted = len(drafted_tokens)
num_accepted = 0

print("--- Simulating Speculative Decoding with Identical Models ---")

for i, token in enumerate(drafted_tokens):
    # For each drafted token, get the probability from both models.
    # Since the models are identical, these probabilities will be the same.
    p_draft = draft_model(context)[token]
    p_target = target_model(context)[token]

    # The acceptance ratio is P_target / P_draft.
    acceptance_ratio = p_target / p_draft

    # The token is accepted if a random number in [0, 1) is less than min(1, acceptance_ratio).
    # Since acceptance_ratio is 1, min(1, 1) is 1. A random number from [0, 1)
    # is always less than 1, so acceptance is guaranteed.
    is_accepted = np.random.random() < min(1, acceptance_ratio)

    # In a real implementation, `is_accepted` would be True due to the math,
    # barring any floating point inaccuracies. We expect it to be true.
    if is_accepted:
        num_accepted += 1

    print(f"\nStep {i+1}: Evaluating token '{token}'")
    print(f"  - Probability from Draft Model:  {p_draft}")
    print(f"  - Probability from Target Model: {p_target}")
    print(f"  - Acceptance Ratio (P_target / P_draft): {p_target} / {p_draft} = {acceptance_ratio}")
    print(f"  - Result: Token Accepted")

# Finally, we calculate the overall acceptance rate.
final_acceptance_rate = num_accepted / num_drafted

print("\n--- Final Calculation ---")
print(f"Total number of accepted tokens: {num_accepted}")
print(f"Total number of drafted tokens:  {num_drafted}")
print(f"Final Acceptance Rate = {num_accepted} / {num_drafted} = {final_acceptance_rate}")
