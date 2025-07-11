def calculate_min_passes():
    """
    Calculates the minimal number of forward passes for a specific
    multiple-choice question scenario.
    """
    # Number of choices of each type
    num_single_token_choices = 4
    num_multi_token_choices = 4

    # For the minimal case, multi-token answers must have the shortest
    # possible length, which is 2 ("more than one").
    min_len_multi_token = 2

    # A single forward pass on the prompt is sufficient to get the likelihoods
    # for all single-token answers. This pass ALSO gives the likelihood of the
    # FIRST token for all multi-token answers.
    initial_shared_pass = 1

    # For each multi-token choice of length k, we need (k-1) additional passes.
    # In our minimal case, k=2, so we need 1 additional pass per choice.
    additional_passes_per_multi_choice = min_len_multi_token - 1
    total_additional_passes = num_multi_token_choices * additional_passes_per_multi_choice

    # The total is the sum of the initial shared pass and all the additional passes.
    total_minimal_passes = initial_shared_pass + total_additional_passes

    print("To find the answer with the highest log-likelihood, we must calculate it for all 8 choices.")
    print("-" * 30)
    print(f"1. One initial forward pass is run on the prompt. This single pass gives us:")
    print(f"   - The complete log-likelihood for all {num_single_token_choices} single-token answers.")
    print(f"   - The log-likelihood of the FIRST token for all {num_multi_token_choices} multi-token answers.")
    print(f"Initial shared passes: {initial_shared_pass}")
    print()

    print(f"2. To find the minimal total passes, we assume each of the {num_multi_token_choices} multi-token answers has the minimum possible length of {min_len_multi_token} tokens.")
    print(f"   - For each of these, we need to find the likelihood of the second token.")
    print(f"   - This requires {additional_passes_per_multi_choice} additional pass for each of the {num_multi_token_choices} choices.")
    print(f"Total additional passes needed: {total_additional_passes}")
    print()

    print("3. The total minimal number of passes is the sum of the shared pass and the additional passes.")
    print("Final Equation:")
    # We explicitly show each additional pass to be clear
    equation_parts = [str(initial_shared_pass)] + [str(additional_passes_per_multi_choice)] * num_multi_token_choices
    print(f"  {' + '.join(equation_parts)} = {total_minimal_passes}")


if __name__ == "__main__":
    calculate_min_passes()