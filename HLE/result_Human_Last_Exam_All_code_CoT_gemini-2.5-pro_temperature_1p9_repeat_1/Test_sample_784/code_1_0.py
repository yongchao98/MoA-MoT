def calculate_forward_passes():
    """
    Calculates the minimal number of forward passes for a specific
    multiple-choice question setup.
    """
    # Number of answer choices of each type
    num_single_token_choices = 4
    num_multi_token_choices = 4

    print("--- Calculating Passes for Single-Token Choices ---")
    # For all choices that are a single token, we can evaluate them simultaneously.
    # A single forward pass on the prompt gives the log likelihoods for all possible
    # next tokens, so we can get the scores for all 4 single-token answers from this one pass.
    passes_for_single_token = 1
    print(f"Number of choices with a single token: {num_single_token_choices}")
    print(f"Forward passes required for all single-token choices: {passes_for_single_token}")
    print("-" * 40)

    print("--- Calculating Passes for Multi-Token Choices ---")
    # For each choice that is a sequence of multiple tokens, we must perform a
    # separate forward pass to calculate its full conditional log likelihood.
    # The context changes for each choice, so they cannot be evaluated in one go.
    passes_for_multi_token = num_multi_token_choices
    print(f"Number of choices with multiple tokens: {num_multi_token_choices}")
    print(f"Forward passes required for each multi-token choice: 1")
    print(f"Total forward passes required for all multi-token choices: {passes_for_multi_token}")
    print("-" * 40)

    # The minimal total is the sum of the two categories.
    total_passes = passes_for_single_token + passes_for_multi_token

    print("--- Total Minimal Forward Passes ---")
    print("The total minimal number of passes is the sum of the passes for each category.")
    print(f"Final Calculation: {passes_for_single_token} + {passes_for_multi_token} = {total_passes}")


if __name__ == '__main__':
    calculate_forward_passes()
    print("<<<5>>>")