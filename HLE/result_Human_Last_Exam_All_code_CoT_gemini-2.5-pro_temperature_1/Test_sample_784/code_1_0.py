def calculate_min_forward_passes():
    """
    Calculates the minimal number of forward passes for an 8-choice MCQ.
    
    The question has:
    - 4 choices with a single output token.
    - 4 choices with more than one output token.
    """
    
    # Number of choices with a single token
    num_single_token_choices = 4

    # Number of choices with more than one token
    num_multi_token_choices = 4

    # To find the MINIMAL number of passes, we assume the minimal possible length
    # for the multi-token choices. Since they must have "more than one" token,
    # the minimal length is 2.
    min_len_multi_token = 2

    # --- Step 1: The Initial Pass ---
    # A single forward pass is run with the question prompt as input.
    # This one pass is sufficient to calculate the full probability for all single-token choices.
    # It also gives us the probability of the *first* token for all multi-token choices.
    initial_pass = 1
    
    print("### Calculation Breakdown ###")
    print(f"1. An initial forward pass is required on the question prompt.")
    print(f"   - This single pass provides the probabilities needed for all {num_single_token_choices} single-token choices.")
    print(f"   - It also provides the probability of the *first* token for the {num_multi_token_choices} multi-token choices.")
    print(f"   - Passes used so far: {initial_pass}\n")

    # --- Step 2: Additional Passes for Multi-Token Choices ---
    # For each multi-token choice of length L, we need L-1 additional passes after the initial one.
    # Since we assume the minimal length is 2, we need 2 - 1 = 1 additional pass per multi-token choice.
    additional_passes_per_choice = min_len_multi_token - 1
    total_additional_passes = num_multi_token_choices * additional_passes_per_choice

    print(f"2. Additional passes are needed for the {num_multi_token_choices} multi-token choices.")
    print(f"   - To ensure the number is *minimal*, we assume each has the shortest possible length: {min_len_multi_token} tokens.")
    print(f"   - After the initial pass, each of these choices requires {additional_passes_per_choice} more pass(es) to evaluate its second token.")
    print(f"   - Total additional passes = {num_multi_token_choices} choices * {additional_passes_per_choice} pass/choice = {total_additional_passes}\n")

    # --- Step 3: Total Calculation ---
    total_passes = initial_pass + total_additional_passes
    
    print("3. The minimal total number of forward passes is the sum.")
    print(f"   Final Equation: {initial_pass} (initial pass) + {total_additional_passes} (additional passes) = {total_passes}")

calculate_min_forward_passes()