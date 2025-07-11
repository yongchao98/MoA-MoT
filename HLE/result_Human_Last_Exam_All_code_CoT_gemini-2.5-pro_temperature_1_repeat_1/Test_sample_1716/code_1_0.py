import sys

# Step 1: Define Alice's payoffs for the symmetric outcomes.
# The problem states that Alice and Bob are superrational agents. This means
# they know that since their reasoning is identical, they will make the same choice.
# We only need to consider the cases where both Alice and Bob do the same thing.

# If both Rest, Alice's payoff is 0.
payoff_both_rest = 0

# If both Bike, Alice's payoff is -2.
payoff_both_bike = -2

# If both Run, Alice's payoff is -3.
payoff_both_run = -3

# Step 2: Explain the reasoning.
# A superrational agent aims to maximize their utility, knowing the other agent
# will do the same, leading to a symmetric outcome.
print("As superrational agents, Alice and Bob know they will make the same choice.")
print("Therefore, we evaluate Alice's utility for each possible symmetric outcome:")
print(f"  - If both Rest, Alice's utility is {payoff_both_rest}.")
print(f"  - If both Bike, Alice's utility is {payoff_both_bike}.")
print(f"  - If both Run, Alice's utility is {payoff_both_run}.")
print("")

# Step 3: Determine the optimal choice.
# Alice will choose the action that leads to the best outcome for herself from
# among these possibilities. This involves a simple comparison.
utilities = [payoff_both_rest, payoff_both_bike, payoff_both_run]
expected_utility = max(utilities)

# Step 4: Print the final calculation and result.
# The 'equation' is the comparison of the possible payoffs.
print("To maximize her utility, Alice chooses the action with the highest payoff.")
# Using str.join to format the numbers for the max function presentation
utility_numbers_str = ", ".join(map(str, utilities))
print(f"Alice's Expected Utility = max({utility_numbers_str})")
print(f"Alice's Expected Utility = {expected_utility}")

# The user has asked for the final answer in a specific format.
# However, the instructions for me say "don't ask users to copy and paste the result".
# I will output the final answer in the requested format as the last line of my response,
# separate from the code block, as per the examples.
# The code itself will print the result as its final step.
# Redirecting the final answer to a different stream to comply with the format <<<answer>>>
# while keeping the python code runnable and self-contained.
sys.stderr.write(f"<<<{expected_utility}>>>")