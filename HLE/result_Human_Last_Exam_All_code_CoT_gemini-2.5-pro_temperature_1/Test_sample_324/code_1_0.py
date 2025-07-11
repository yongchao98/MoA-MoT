import torch
import torch.nn.functional as F

# This is a simplified simulation to demonstrate the core logic.
# In this sanity check, the draft and target models are the same.
# We can represent them with a single function that returns a probability distribution.

def model_logic(input_sequence):
    """
    A mock, deterministic model function. Given an input sequence,
    it returns a fixed probability distribution (logits) for the next token.
    Using a seed based on the input ensures it's deterministic.
    """
    torch.manual_seed(sum(input_sequence) + len(input_sequence))
    # In a real LLM, this would be a complex transformer forward pass.
    # The output dimension represents the vocabulary size.
    vocab_size = 50 
    logits = torch.randn(vocab_size)
    return logits

def get_probabilities(logits):
    """Convert logits to probabilities using softmax."""
    return F.softmax(logits, dim=-1)

# --- Simulation ---

# Configuration
initial_prompt = [10, 25, 3]
num_draft_tokens = 4 # The number of tokens the draft model will generate

# Variables to track the results
draft_tokens = []
accepted_count = 0
current_sequence = list(initial_prompt)

# 1. DRAFTING PHASE: The "draft model" generates a sequence of tokens.
print(f"--- Sanity Check: Draft Model == Target Model ---")
print(f"Initial Prompt: {initial_prompt}")
print(f"Drafting {num_draft_tokens} tokens...\n")
temp_sequence = list(current_sequence)
for _ in range(num_draft_tokens):
    # Get logits from the draft model
    draft_logits = model_logic(temp_sequence)
    # For demonstration, we'll greedily pick the most likely token
    next_token = torch.argmax(draft_logits).item()
    draft_tokens.append(next_token)
    temp_sequence.append(next_token)

print(f"Drafted Tokens: {draft_tokens}\n")
print(f"--- Verification Phase ---")

# 2. VERIFICATION PHASE: We check each drafted token against the "target model".
for i, token in enumerate(draft_tokens):
    # The context for verifying the current token
    context = current_sequence + draft_tokens[:i]
    
    # Since draft_model == target_model, they use the same logic
    draft_logits = model_logic(context)
    target_logits = model_logic(context)
    
    # The probabilities will therefore be identical
    draft_probs = get_probabilities(draft_logits)
    target_probs = get_probabilities(target_logits)
    
    prob_from_draft_model = draft_probs[token]
    prob_from_target_model = target_probs[token]
    
    # Acceptance probability is min(1, P_target / P_draft)
    # Since P_target == P_draft, the ratio is 1.
    acceptance_prob = min(1.0, prob_from_target_model / prob_from_draft_model)
    
    print(f"Verifying token '{token}' at position {i+1}:")
    print(f"  P_target({token}) / P_draft({token}) = {prob_from_target_model:.4f} / {prob_from_draft_model:.4f} = 1.0")
    print(f"  Acceptance probability = min(1, 1.0) = {acceptance_prob:.1f}. ACCEPTED.")
    
    if acceptance_prob == 1.0:
        accepted_count += 1
    else:
        # This part will not be reached in this scenario
        break

# 3. FINAL CALCULATION
total_drafted = len(draft_tokens)
acceptance_rate = accepted_count / total_drafted

print("\n--- Final Result ---")
print(f"Total tokens drafted: {total_drafted}")
print(f"Total tokens accepted: {accepted_count}")
print(f"Final acceptance rate equation: {accepted_count} / {total_drafted} = {acceptance_rate}")