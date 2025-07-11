# Define the number of choices of each type.
num_multi_token_choices = 4

# A single forward pass is needed to get the probability distribution for the first token.
# This one pass covers the 4 single-token answers and the first token of the 4 multi-token answers.
passes_for_first_token = 1

# To find the minimal number of passes, we assume the multi-token answers have the minimum
# possible number of tokens, which is 2. We've already accounted for the first token.
# We need one additional forward pass for each of the 4 multi-token answers to evaluate their second token.
additional_passes_for_multi_token_choices = num_multi_token_choices

# The total minimal number of forward passes is the sum.
total_passes = passes_for_first_token + additional_passes_for_multi_token_choices

# Print the final equation showing how the total is calculated.
print(f"{passes_for_first_token} + {additional_passes_for_multi_token_choices} = {total_passes}")