def calculate_min_forward_passes():
    """
    Calculates the minimal number of forward passes for an 8-choice MCQ.

    The problem states:
    - 8 total choices
    - 4 choices are single-token
    - 4 choices are multi-token
    """
    
    single_token_choices = 4
    multi_token_choices = 4
    
    # Step 1: Passes for single-token choices
    # To get the log likelihood for any single token, we run one forward pass with the
    # prompt. The output is a probability distribution over the entire vocabulary for the
    # next token. We can look up the probabilities for all 4 single-token choices
    # from this one result.
    passes_for_single_token = 1
    
    print("Plan to calculate the minimal number of forward passes:")
    print("1. Calculate passes for the 4 single-token choices.")
    print("2. Calculate passes for the 4 multi-token choices.")
    print("3. Sum the results for the total.\n")
    
    print(f"Step 1: Evaluating the {single_token_choices} single-token choices.")
    print("A single forward pass on the prompt is sufficient to get the log likelihoods for all possible next tokens.")
    print(f"Passes needed for this step: {passes_for_single_token}\n")
    
    # Step 2: Passes for multi-token choices
    # Each multi-token choice is a unique sequence. To calculate the conditional
    # log likelihood of a full sequence, it must be fed to the model in its own
    # forward pass. Therefore, each of the 4 multi-token choices requires
    # its own pass.
    passes_for_multi_token = multi_token_choices
    
    print(f"Step 2: Evaluating the {multi_token_choices} multi-token choices.")
    print("Each multi-token sequence requires a separate forward pass to calculate its full log likelihood.")
    print(f"Passes needed for this step: {passes_for_multi_token}\n")

    # Step 3: Total passes
    # The total minimal number of passes is the sum from both steps.
    total_passes = passes_for_single_token + passes_for_multi_token
    
    print("Step 3: Calculating the total.")
    print("The total minimal number of forward passes is the sum of the passes from each step.")
    print(f"Total passes = (passes for single-token choices) + (passes for multi-token choices)")
    print(f"Final Equation: {passes_for_single_token} + {passes_for_multi_token} = {total_passes}")

calculate_min_forward_passes()
<<<5>>>