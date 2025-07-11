def calculate_minimal_forward_passes():
    """
    Calculates the minimal number of forward passes for an 8-choice MCQ
    with a specific mix of single- and multi-token answers.
    """
    # There are 4 answer choices that consist of a single output token.
    num_single_token_choices = 4

    # There are 4 answer choices that consist of more than one output token.
    num_multi_token_choices = 4

    # For all single-token choices, we can run one forward pass on the prompt
    # and find their respective probabilities from the resulting next-token
    # probability distribution.
    passes_for_single_token_choices = 1

    # For each multi-token choice, we must run a separate forward pass on the
    # sequence of [prompt + answer_choice] to calculate its full probability.
    # Therefore, the number of passes equals the number of choices.
    passes_for_multi_token_choices = num_multi_token_choices

    # The total minimal number of passes is the sum of the two cases.
    total_passes = passes_for_single_token_choices + passes_for_multi_token_choices

    print(f"Number of single-token choices: {num_single_token_choices}")
    print(f"Number of multi-token choices: {num_multi_token_choices}")
    print("-" * 30)
    print(f"Passes needed for all single-token choices: {passes_for_single_token_choices}")
    print(f"Passes needed for all multi-token choices: {passes_for_multi_token_choices}")
    print("-" * 30)
    print("The minimal total number of forward passes is the sum.")
    # Final equation showing each number
    print(f"Total Passes = {passes_for_single_token_choices} + {passes_for_multi_token_choices} = {total_passes}")

calculate_minimal_forward_passes()