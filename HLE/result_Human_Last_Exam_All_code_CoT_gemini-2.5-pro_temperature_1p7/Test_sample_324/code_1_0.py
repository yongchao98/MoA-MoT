import numpy as np

def simulated_llm(context_tokens):
    """
    This function simulates an LLM. For a given context, it deterministically
    produces a probability distribution for the next token. Since the draft and
    target models are the same in this problem, we only need one such function.
    """
    # Use the context to seed the random number generator for deterministic output
    # based on input. The hash ensures we get a valid integer seed.
    seed = hash(tuple(context_tokens)) & 0xFFFFFFFF
    np.random.seed(seed)
    
    # Generate a dummy probability distribution over a vocabulary of 500 tokens
    vocab_size = 500
    probabilities = np.random.rand(vocab_size)
    return probabilities / np.sum(probabilities) # Normalize to sum to 1

def run_speculative_decoding_sanity_check(draft_length_k=5):
    """
    Simulates speculative decoding where the draft model and target model are identical.
    """
    # 1. Setup
    initial_context = [15, 83, 101] # An example initial prompt
    
    # In this scenario, both models are the same function
    draft_model = simulated_llm
    target_model = simulated_llm

    # 2. Draft Phase: The draft model generates 'k' tokens
    drafted_tokens = []
    current_context = list(initial_context)
    for _ in range(draft_length_k):
        # Get probability distribution from the draft model
        prob_dist = draft_model(current_context)
        # Choose the most likely token (greedy decoding)
        next_token = np.argmax(prob_dist)
        drafted_tokens.append(next_token)
        current_context.append(next_token)

    # 3. Verification Phase: The target model verifies the drafted tokens
    accepted_count = 0
    verification_context = list(initial_context)
    for i in range(draft_length_k):
        drafted_token = drafted_tokens[i]

        # Get the probability distributions from both models for the current step.
        # Since the models and context are identical, these distributions will be identical.
        p_dist_draft = draft_model(verification_context)
        p_dist_target = target_model(verification_context)
        
        # Get the specific probability of the token that was drafted
        p_draft = p_dist_draft[drafted_token]
        p_target = p_dist_target[drafted_token]

        # Since p_draft and p_target are identical, their ratio is 1.0
        acceptance_prob = min(1.0, p_target / p_draft)

        # The acceptance check will always pass
        if np.random.random() < acceptance_prob:
            accepted_count += 1
            verification_context.append(drafted_token)
        else:
            # This 'else' block will never be reached in this scenario
            break
            
    # 4. Results
    total_drafted_tokens = len(drafted_tokens)
    final_acceptance_rate = accepted_count / total_drafted_tokens

    print(f"Speculative Decoding Sanity Check Results:")
    print(f"Total tokens drafted: {total_drafted_tokens}")
    print(f"Total tokens accepted: {accepted_count}")
    print("\nFinal acceptance rate equation:")
    print(f"{accepted_count} / {total_drafted_tokens} = {final_acceptance_rate}")
    print("\nAs expected, the acceptance rate is precisely 1 because the draft and target models are identical.")


if __name__ == '__main__':
    # Run the simulation
    run_speculative_decoding_sanity_check()