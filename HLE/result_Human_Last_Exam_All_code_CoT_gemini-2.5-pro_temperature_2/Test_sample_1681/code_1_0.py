import numpy as np

# Payoff matrix values
R = 5  # Reward for mutual cooperation (Cooperate, Cooperate)
S = 0  # Sucker's payoff (Cooperate, Defect)
T = 7  # Temptation to defect (Defect, Cooperate)
P = 1  # Punishment for mutual defection (Defect, Defect)

# In the EWL quantum game model with optimal entanglement and full SU(2) strategies,
# a new Nash Equilibrium emerges where both players pick the "Quantum" strategy Q.
# Our analysis shows that playing (Q, Q) results in the same outcome as mutual
# cooperation (C, C) with 100% certainty.

# Therefore, the probabilities of the four classical outcomes are:
P_CC = 1.0  # Probability of (Cooperate, Cooperate)
P_CD = 0.0  # Probability of (Cooperate, Defect)
P_DC = 0.0  # Probability of (Defect, Cooperate)
P_DD = 0.0  # Probability of (Defect, Defect)

# The payoff for each player at this equilibrium is calculated from the
# expected value based on these probabilities and the payoff matrix.
equilibrium_payoff = P_CC * R + P_CD * S + P_DC * T + P_DD * P

# Display the final equation and the result
print("The equilibrium payoff is found by calculating the expected value for a player at the (Q, Q) Nash Equilibrium.")
print("The outcome probabilities at this equilibrium are P_CC=1.0, P_CD=0.0, P_DC=0.0, P_DD=0.0.")
print("\nThe final equation for a player's payoff is:")
print(f"Payoff = P_CC * R + P_CD * S + P_DC * T + P_DD * P")
print(f"Payoff = {P_CC} * {R} + {P_CD} * {S} + {P_DC} * {T} + {P_DD} * {P} = {equilibrium_payoff}")

print(f"\nThe quantum equilibrium point is a payoff of {equilibrium_payoff} for each player, resolving the dilemma.")
