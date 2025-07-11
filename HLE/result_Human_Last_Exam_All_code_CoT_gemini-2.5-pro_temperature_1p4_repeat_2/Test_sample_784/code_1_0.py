def calculate_min_forward_passes():
    """
    Calculates the minimal number of forward passes for an 8-choice MCQ.
    The problem states 4 choices are single-token and 4 are multi-token.
    """

    # --- Problem Definition ---
    num_single_token_choices = 4
    num_multi_token_choices = 4

    print("To find the minimal number of forward passes, we must assume the most optimal scenario for evaluation.\n")

    # --- Step 1: The Initial Forward Pass ---
    # A single forward pass on the main question prompt is always required.
    # This pass provides the probabilities for the very next token.
    # This is sufficient to:
    #   a) Evaluate all 4 single-token choices.
    #   b) Evaluate the *first* token of all 4 multi-token choices.
    passes_for_prompt = 1
    print(f"Pass 1: A forward pass on the question prompt.")
    print(f"- This single pass evaluates all {num_single_token_choices} single-token choices.")
    print(f"- It also gives the probability for the first token of all {num_multi_token_choices} multi-token choices.")
    print(f"Passes used so far: {passes_for_prompt}\n")

    # --- Step 2: Passes for Multi-Token Choices (Optimal Case) ---
    # To minimize passes, we assume the best case for the multi-token choices:
    #  a) They are all the minimum possible length for a multi-token choice (2 tokens).
    #  b) They all share the same first token.
    #
    # Because they all share the same first token, we only need one additional
    # forward pass on the shared prefix (prompt + shared_first_token) to get
    # the probabilities for the second token of all 4 choices.
    additional_passes_for_multi_token = 1
    print(f"Pass 2: A forward pass on the shared multi-token prefix.")
    print(f"In the optimal case, all {num_multi_token_choices} multi-token choices share the same first token.")
    print(f"Therefore, only {additional_passes_for_multi_token} more pass is needed to evaluate the second token of all of them.")
    print(f"Additional passes needed: {additional_passes_for_multi_token}\n")

    # --- Step 3: Total Calculation ---
    total_passes = passes_for_prompt + additional_passes_for_multi_token
    print("The total minimal number of passes is the sum of the initial pass and the additional pass for the shared multi-token prefix.")
    print(f"Final Equation: {passes_for_prompt} (for prompt) + {additional_passes_for_multi_token} (for shared prefix)")
    print(f"Result: {passes_for_prompt} + {additional_passes_for_multi_token} = {total_passes}")

    return total_passes

if __name__ == "__main__":
    calculate_min_forward_passes()