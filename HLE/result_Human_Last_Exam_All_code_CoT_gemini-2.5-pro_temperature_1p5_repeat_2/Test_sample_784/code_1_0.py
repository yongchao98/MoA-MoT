def solve_forward_pass_query():
    """
    Calculates the minimal number of forward passes for an 8-choice MCQ.

    The problem states there are 8 choices in total:
    - 4 choices consist of a single output token.
    - 4 choices consist of more than one output token.

    The goal is to find the chosen answer, defined as the one with the
    highest conditional log likelihood, and to do so with the minimum
    number of forward passes.
    """

    # Number of choices of each type
    num_single_token_choices = 4
    num_multi_token_choices = 4

    # Step 1: Calculate passes for single-token choices.
    # To get the log likelihoods for all single-token answers (e.g., log P(choice | question)),
    # we can run a single forward pass with the question as input. The output provides the
    # probability distribution for the very next token, from which we can find the
    # likelihoods for all 4 single-token options.
    passes_for_single_token_choices = 1

    # Step 2: Calculate passes for multi-token choices.
    # For each multi-token answer, we need to calculate the likelihood of the entire
    # sequence. For a sequence (t1, t2, ..., tk), this is calculated as the sum of
    # log P(t1|q), log P(t2|q,t1), etc. A causal transformer calculates this for
    # one specific sequence in a single forward pass. Since each of the 4 multi-token
    # answers is a different sequence, each requires its own forward pass.
    passes_for_multi_token_choices = num_multi_token_choices

    # Step 3: Sum the passes for the final result.
    total_passes = passes_for_single_token_choices + passes_for_multi_token_choices

    # Print the detailed breakdown and the final calculation.
    print("Calculating the minimal number of forward passes for an 8-choice MCQ:")
    print(f" - {num_single_token_choices} choices are single-token.")
    print(f" - {num_multi_token_choices} choices are multi-token.")
    print("\nBreakdown of passes:")
    print(f"1. A single forward pass is needed to evaluate all {num_single_token_choices} single-token choices simultaneously.")
    print(f"2. A separate forward pass is needed for each of the {num_multi_token_choices} multi-token choices.")
    print("\nThe final calculation is:")
    print(f"Total Passes = {passes_for_single_token_choices} (for single-token) + {passes_for_multi_token_choices} (for multi-token) = {total_passes}")

solve_forward_pass_query()