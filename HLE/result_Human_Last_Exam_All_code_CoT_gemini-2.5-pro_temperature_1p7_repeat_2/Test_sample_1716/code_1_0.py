# Step 1: Define the payoff matrix for Alice based on the problem description.
# The payoffs are structured as a dictionary where the first key is Alice's move
# and the second key is Bob's move.
# Payoffs are listed as (Alice, Bob) in the problem, we only need Alice's here.
alice_payoffs = {
    # Alice's Move: Bob's Move -> Alice's Payoff
    'Rest': {'Rest': 0, 'Bike': 2, 'Run': 4},
    'Bike': {'Rest': 0, 'Bike': -2, 'Run': 2},
    'Run':  {'Rest': 0, 'Bike': 0, 'Run': -3}
}

# Step 2 & 3: Apply the concept of superrationality.
# Alice and Bob, as superrational agents in a symmetric game, can coordinate on the
# optimal symmetric correlated equilibrium. This allows them to achieve a higher
# joint payoff than simply picking the same action.
# The optimal equilibrium in this case is to randomly choose between two scenarios
# with equal probability, ensuring fairness and optimality:
# Scenario A: Alice Rests, Bob Runs
# Scenario B: Alice Runs, Bob Rests

# The probability for each of these scenarios is 0.5.
prob_scenario_A = 0.5
prob_scenario_B = 0.5

# Step 4: Calculate Alice's expected utility.
# Get Alice's payoff from Scenario A (Alice Rests, Bob Runs).
payoff_A = alice_payoffs['Rest']['Run']

# Get Alice's payoff from Scenario B (Alice Runs, Bob Rests).
payoff_B = alice_payoffs['Run']['Rest']

# The expected utility is the sum of each payoff multiplied by its probability.
expected_utility = (prob_scenario_A * payoff_A) + (prob_scenario_B * payoff_B)

print("To find Alice's expected utility, we consider the optimal strategy for superrational agents.")
print("This is a symmetric correlated equilibrium with two outcomes, each occurring with 50% probability:\n")
print(f"1. Alice Rests and Bob Runs: Alice's payoff is {payoff_A}")
print(f"2. Alice Runs and Bob Rests: Alice's payoff is {payoff_B}\n")
print("Alice's expected utility is calculated as follows:")
print(f"({prob_scenario_A} * {payoff_A}) + ({prob_scenario_B} * {payoff_B}) = {expected_utility}")
