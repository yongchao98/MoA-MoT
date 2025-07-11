def calculate_min_forward_passes():
    """
    Calculates the minimum number of forward passes for an 8-choice MCQ.
    The MCQ has 4 single-token and 4 multi-token answer choices.
    """
    # Define the parameters of the problem
    num_single_token_choices = 4
    num_multi_token_choices = 4
    
    # To find the 'minimal' number of passes, we assume the multi-token choices
    # have the minimum possible number of tokens greater than one, which is 2.
    min_tokens_in_multi_token_choice = 2

    # --- Step 1: The Initial Forward Pass ---
    # A single forward pass is run on the input prompt. The output of this pass
    # is a probability distribution over the entire vocabulary for the next token.
    # This single pass is sufficient to determine the probability for all 4 single-token
    # choices, as well as the *first* token of all 4 multi-token choices.
    initial_pass = 1
    
    # --- Step 2: Additional Passes for Multi-Token Choices ---
    # For a multi-token answer with 'k' tokens (t_1, t_2, ..., t_k), the probability is:
    # P(t_1|prompt) * P(t_2|prompt,t_1) * ... * P(t_k|prompt,t_1,...,t_{k-1})
    # Calculating P(t_1|prompt) is done in the initial pass.
    # Each subsequent token (t_2, t_3, etc.) requires its own new forward pass
    # with the preceding tokens added to the input.
    # Therefore, a k-token answer requires k-1 additional passes.
    additional_passes_per_multi_choice = min_tokens_in_multi_token_choice - 1
    
    # --- Step 3: Total Calculation ---
    # The total number of passes is the initial pass plus the sum of all
    # additional passes required for each of the multi-token choices.
    total_additional_passes = num_multi_token_choices * additional_passes_per_multi_choice
    
    total_passes = initial_pass + total_additional_passes
    
    # --- Step 4: Output the explanation and result ---
    print("Calculating the minimal number of forward passes:")
    print("-" * 50)
    print(f"An initial forward pass is required for the prompt. This covers all {num_single_token_choices} single-token choices and the first token of all {num_multi_token_choices} multi-token choices.")
    print(f"Initial pass = {initial_pass}")
    print("-" * 50)
    print(f"For each of the {num_multi_token_choices} multi-token choices, we assume the minimal length of {min_tokens_in_multi_token_choice} tokens.")
    print(f"Each multi-token choice needs ({min_tokens_in_multi_token_choice} - 1) = {additional_passes_per_multi_choice} additional pass.")
    print(f"Total additional passes = {num_multi_token_choices} choices * {additional_passes_per_multi_choice} pass/choice = {total_additional_passes}")
    print("-" * 50)
    print("The final calculation for the total minimal number of passes is:")
    print(f"Total Passes = (Initial Pass) + (Passes for 2nd token of choice 5) + (Passes for 2nd token of choice 6) + (Passes for 2nd token of choice 7) + (Passes for 2nd token of choice 8)")
    print(f"Total Passes = {initial_pass} + {additional_passes_per_multi_choice} + {additional_passes_per_multi_choice} + {additional_passes_per_multi_choice} + {additional_passes_per_multi_choice}")
    print(f"Total Passes = {total_passes}")

calculate_min_forward_passes()