# Define Alice's payoffs for the symmetric outcomes where both agents choose the same action.
# If both rest, the payoff is 0.
payoff_rest = 0
# If both bike, the payoff is -2.
payoff_bike = -2
# If both run, the payoff is -3.
payoff_run = -3

# Superrational agents know they will both choose the same action.
# Therefore, they choose the action that gives the best payoff on the diagonal of the payoff matrix.
# We find the maximum of these potential outcomes.
alices_utility = max(payoff_rest, payoff_bike, payoff_run)

# Print the explanation and the final equation for Alice's expected utility.
print("Because Alice and Bob are superrational, they know they will both make the same choice.")
print("Alice's decision is based on choosing the action that provides the best outcome if both players choose it.")
print("The possible payoffs for Alice are:")
print(f"  - Both Rest: {payoff_rest}")
print(f"  - Both Bike: {payoff_bike}")
print(f"  - Both Run: {payoff_run}")
print("\nTo find her expected utility, Alice calculates the maximum of these values:")
print(f"Alice's Expected Utility = max({payoff_rest}, {payoff_bike}, {payoff_run}) = {alices_utility}")