import numpy as np

# Step 1: Define the choices and construct the payoff matrix from the story.
# The choices are Rest, Bike, and Run.
choices = ["Rest", "Bike", "Run"]

# The payoff matrix is structured as (Alice's payoff, Bob's payoff).
# Rows represent Alice's choice, Columns represent Bob's choice.
payoffs = np.array([
    # Bob chooses: Rest      Bike        Run
    [(0, 0),    (2, 0),     (4, 0)],   # Alice chooses Rest
    [(0, 2),    (-2, -2),   (2, 0)],   # Alice chooses Bike
    [(0, 4),    (0, 2),     (-3, -3)]  # Alice chooses Run
])

# We can separate the payoffs for clarity.
alice_payoffs = payoffs[:, :, 0]
bob_payoffs = payoffs[:, :, 1]

# Step 2 & 3: Analyze the game by finding weakly dominated strategies.

# --- Analysis for Alice ---
print("Analyzing for Alice's strategies:")
alice_dominated_indices = set()
# For each strategy `i`, we check if any other strategy `j` dominates it.
for i in range(len(choices)):
    for j in range(len(choices)):
        if i == j:
            continue
        
        # Strategy `j` weakly dominates `i` if:
        # 1. Payoff(j) is never worse than Payoff(i) for any of Bob's moves.
        # 2. Payoff(j) is strictly better than Payoff(i) for at least one of Bob's moves.
        is_never_worse = all(alice_payoffs[j, k] >= alice_payoffs[i, k] for k in range(len(choices)))
        is_sometimes_better = any(alice_payoffs[j, k] > alice_payoffs[i, k] for k in range(len(choices)))
        
        if is_never_worse and is_sometimes_better:
            print(f"- Alice's choice '{choices[j]}' weakly dominates '{choices[i]}'.")
            alice_dominated_indices.add(i)

# --- Analysis for Bob ---
print("\nAnalyzing for Bob's strategies:")
bob_dominated_indices = set()
# Bob is the column player. We compare columns `j` and `j_prime`.
for j in range(len(choices)):
    for j_prime in range(len(choices)):
        if j == j_prime:
            continue

        # Strategy `j_prime` weakly dominates `j` if:
        # 1. Payoff(j_prime) is never worse than Payoff(j) for any of Alice's moves.
        # 2. Payoff(j_prime) is strictly better than Payoff(j) for at least one of Alice's moves.
        is_never_worse = all(bob_payoffs[k, j_prime] >= bob_payoffs[k, j] for k in range(len(choices)))
        is_sometimes_better = any(bob_payoffs[k, j_prime] > bob_payoffs[k, j] for k in range(len(choices)))
        
        if is_never_worse and is_sometimes_better:
            print(f"- Bob's choice '{choices[j_prime]}' weakly dominates '{choices[j]}'.")
            bob_dominated_indices.add(j)
            
# Step 4: Deduce the outcome based on superrationality.
# Superrational players will not use weakly dominated strategies.
alice_undominated_indices = [i for i in range(len(choices)) if i not in alice_dominated_indices]
bob_undominated_indices = [j for j in range(len(choices)) if j not in bob_dominated_indices]

# Since there is only one undominated strategy for each, they will choose it.
alice_final_choice_index = alice_undominated_indices[0]
bob_final_choice_index = bob_undominated_indices[0]

alice_final_choice = choices[alice_final_choice_index]
bob_final_choice = choices[bob_final_choice_index]

print(f"\nRational deduction leads to a single choice for each player:")
print(f"Alice's only undominated choice is '{alice_final_choice}'.")
print(f"Bob's only undominated choice is '{bob_final_choice}'.")

# Step 5: Calculate Alice's expected utility from the determined outcome.
final_payoffs = payoffs[alice_final_choice_index, bob_final_choice_index]
alice_final_utility = final_payoffs[0]

print("\nThe final outcome is (Alice: {}, Bob: {}).".format(alice_final_choice, bob_final_choice))
print("The resulting payoffs are (Alice: {}, Bob: {}).".format(final_payoffs[0], final_payoffs[1]))

print("\nTherefore, Alice's expected utility is determined by this outcome.")
print("Final Utility Equation:")
print(f"Alice's Utility = Payoff_Alice({alice_final_choice}, {bob_final_choice}) = {alice_final_utility}")