def calculate_minimal_forward_passes():
    """
    Calculates the minimal number of forward passes for an 8-choice MCQ.
    
    The breakdown is as follows:
    - 4 choices are single-token.
    - 4 choices are multi-token.
    """
    
    # 1. A single forward pass is needed to evaluate the first token of ALL answer choices.
    #    This single pass is sufficient to get the full log-likelihood for the 4 single-token answers.
    initial_pass_for_first_token = 1
    
    # 2. For each of the 4 multi-token answer choices, we need an additional forward pass for each subsequent token.
    #    To find the MINIMAL number of passes, we assume the multi-token choices are as short as possible,
    #    which is 2 tokens. Therefore, each of these 4 choices requires 1 additional pass for its second token.
    num_multi_token_choices = 4
    additional_passes_per_multi_token_choice = 1 # Assuming 2 tokens total per choice for minimality.
    
    total_additional_passes = num_multi_token_choices * additional_passes_per_multi_token_choice
    
    # 3. The total is the initial pass plus the additional passes for the multi-token choices.
    total_passes = initial_pass_for_first_token + total_additional_passes
    
    print("To find the answer with the highest conditional log-likelihood, we need:")
    print(f"- {initial_pass_for_first_token} forward pass to evaluate all 4 single-token choices and the first token of all 4 multi-token choices.")
    print(f"- {total_additional_passes} additional forward passes (1 for each of the {num_multi_token_choices} multi-token choices to evaluate their second token).")
    print("\nFinal equation:")
    print(f"{initial_pass_for_first_token} + {num_multi_token_choices} = {total_passes}")

calculate_minimal_forward_passes()