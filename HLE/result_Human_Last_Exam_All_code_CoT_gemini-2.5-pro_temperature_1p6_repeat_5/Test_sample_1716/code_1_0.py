# Step 1: Define the payoffs for Alice in the symmetric outcomes.
# These are the only outcomes possible because Alice and Bob are superrational
# and will independently come to the same conclusion.

# If both Alice and Bob rest, her payoff is 0.
utility_both_rest = 0
# If both Alice and Bob bike, her payoff is -2.
utility_both_bike = -2
# If both Alice and Bob run, her payoff is -3.
utility_both_run = -3

# Step 2: As a rational agent, Alice will choose the action that leads
# to the outcome with the highest utility for herself.
alices_utility = max(utility_both_rest, utility_both_bike, utility_both_run)

# Step 3: Print the reasoning and the final calculation.
print("Alice and Bob are superrational, so they will both make the same choice.")
print("Alice evaluates her utility in the three possible symmetric scenarios:")
print(f"1. Both rest: Utility = {utility_both_rest}")
print(f"2. Both bike: Utility = {utility_both_bike}")
print(f"3. Both run: Utility = {utility_both_run}")
print("\nTo maximize her utility, she chooses the best of these outcomes.")
print("The final calculation for her expected utility is:")
print(f"max({utility_both_rest}, {utility_both_bike}, {utility_both_run}) = {alices_utility}")
