# Payoff values from the matrix
# R is the Reward for mutual cooperation (C,C)
# S is the Sucker's payoff for cooperating while the other defects (C,D)
# T is the Temptation to defect while the other cooperates (D,C)
# P is the Punishment for mutual defection (D,D)
R = 5
S = 0
T = 7
P = 1

# In the quantum version with maximal entanglement, a new Nash Equilibrium (Q,Q) emerges.
# When both players play the quantum strategy Q, the probability of the (Cooperate, Cooperate)
# outcome becomes 1, and all other outcomes have a probability of 0.
P_CC = 1.0
P_CD = 0.0
P_DC = 0.0
P_DD = 0.0

# Calculate the payoffs for each player (Alice and Bob)
# Alice's payoff: R*P_CC + S*P_CD + T*P_DC + P*P_DD
# Bob's payoff:   R*P_CC + T*P_CD + S*P_DC + P*P_DD
alice_payoff = R * P_CC + S * P_CD + T * P_DC + P * P_DD
bob_payoff = R * P_CC + T * P_CD + S * P_DC + P * P_DD

# The equilibrium point is the pair of payoffs (Alice's Payoff, Bob's Payoff)
print("The quantum equilibrium resolves the dilemma, leading to the optimal outcome.")
print("The calculation for the equilibrium point is as follows:")

# Print the equation for Alice's payoff
print(f"Alice's Payoff = ({R} * {P_CC}) + ({S} * {P_CD}) + ({T} * {P_DC}) + ({P} * {P_DD}) = {alice_payoff}")

# Print the equation for Bob's payoff
print(f"Bob's Payoff   = ({R} * {P_CC}) + ({T} * {P_CD}) + ({S} * {P_DC}) + ({P} * {P_DD}) = {bob_payoff}")

print(f"\nThe equilibrium point is ({int(alice_payoff)}, {int(bob_payoff)}).")
