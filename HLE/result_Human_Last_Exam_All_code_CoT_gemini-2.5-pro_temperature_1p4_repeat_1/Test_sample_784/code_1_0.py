def calculate_min_forward_passes():
    """
    Calculates the minimal number of forward passes for an 8-choice MCQ.

    The MCQ has:
    - 4 single-token answer choices.
    - 4 multi-token answer choices.
    """

    # Number of choices for each type
    num_single_token_choices = 4
    num_multi_token_choices = 4

    # To find the minimal number of passes, we assume the minimal length
    # for the multi-token answers. "More than one" token means the minimum length is 2.
    min_len_multi_token = 2

    # Step 1: Initial Pass
    # A single forward pass on the prompt gets the probability for the next token.
    # This covers all single-token choices and the first token of all multi-token choices.
    initial_pass = 1

    # Step 2: Additional Passes
    # Each multi-token choice requires additional passes for its subsequent tokens.
    # For a choice of length 2, it needs (2 - 1) = 1 additional pass.
    additional_passes_per_choice = min_len_multi_token - 1

    # In the general case, we need to run these additional passes for each of the
    # multi-token choices separately.
    total_additional_passes = num_multi_token_choices * additional_passes_per_choice

    # Step 3: Total calculation
    total_passes = initial_pass + total_additional_passes

    # --- Outputting the explanation and equation ---
    print("The minimal number of forward passes is calculated as follows:")
    
    # Print the equation part by part to show each number.
    # We start with the initial pass that covers the first token of all choices.
    print(f"\n{initial_pass} (for the first token of all 8 choices)")

    # Add the passes for the remaining tokens of each multi-token choice.
    for i in range(num_multi_token_choices):
        print(f"+ {additional_passes_per_choice} (for the 2nd token of multi-token choice #{i+1})")
    
    # Print the final result
    print(f"\n= {total_passes}")

calculate_min_forward_passes()