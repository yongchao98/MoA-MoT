import random

def calculate_acceptance_rate_for_identical_models():
    """
    Simulates a speculative decoding step where the draft and target models are identical,
    and calculates the expected acceptance rate.
    """

    # K is the number of tokens the draft model proposes in one sequence (the "lookahead").
    K = 5
    
    accepted_tokens_count = 0

    print("--- Sanity Check Simulation: Draft Model == Target Model ---")
    print(f"Simulating one step with K={K} drafted tokens.\n")

    # In a real scenario, the draft model generates K tokens based on a probability distribution.
    # The target model then validates them. We simulate this for each of the K tokens.
    for i in range(K):
        # Let's imagine the draft model picked a token. We don't need to know which one,
        # only the probability it assigned to that token. Let's simulate this probability.
        # This would be P_draft(token | prompt + previous_tokens).
        p_draft = random.uniform(0.01, 1.0)
        
        # Since the draft and target models are identical, the target model will calculate
        # the exact same probability for that same token at that same position.
        # This is the key insight: P_target = P_draft.
        p_target = p_draft

        # The acceptance rule compares the ratio of probabilities.
        # The probability of accepting this token is min(1, P_target / P_draft).
        acceptance_probability = min(1, p_target / p_draft)
        
        # A random number r is drawn from [0, 1). If r < acceptance_probability, the token is accepted.
        # Since acceptance_probability is 1, any r drawn will be less than it. The token is always accepted.
        is_accepted = random.random() < acceptance_probability

        if is_accepted:
            accepted_tokens_count += 1
        
        print(f"Verifying drafted token {i+1}:")
        print(f"  P_draft(token)  = {p_draft:.4f}")
        print(f"  P_target(token) = {p_target:.4f}")
        print(f"  Acceptance Chance = min(1, {p_target:.4f}/{p_draft:.4f}) = {acceptance_probability:.1f}")

    print("\n--- Final Calculation ---")
    
    # The acceptance rate is the ratio of accepted tokens to total drafted tokens.
    final_acceptance_rate = accepted_tokens_count / K

    # In the final output, we need to show the final equation with the numbers.
    print("Final Equation: Acceptance Rate = (Total Accepted Tokens) / (Total Drafted Tokens)")
    print(f"Result: {accepted_tokens_count} / {K} = {final_acceptance_rate}")

# Run the simulation
calculate_acceptance_rate_for_identical_models()