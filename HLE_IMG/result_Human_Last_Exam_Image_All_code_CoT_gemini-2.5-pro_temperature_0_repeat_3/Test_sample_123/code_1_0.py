import numpy as np
from fractions import Fraction

# Let P_X be the probability of reaching the Treasure Room from room X.
# The problem can be modeled by the following system of linear equations,
# where each equation represents the probability of success from a given room
# based on the probabilities of success from the rooms it connects to.
print("The system of equations for the probabilities (P_X) is:")
print("P_ME = (1/2) * P_GH + (1/2) * P_CK")
print("P_GH = (1/4) * P_ME + (1/4) * P_SP + (1/4) * P_L + (1/4) * P_KH")
print("P_CK = (1/2) * P_ME + (1/2) * P_KH")
print("P_SP = (1/4) * P_GH + (1/4) * P_TC + (1/4) * 1 + (1/4) * 0")
print("P_L  = (1/3) * P_GH + (1/3) * P_KH + (1/3) * P_QR")
print("P_KH = (1/4) * P_GH + (1/4) * P_CK + (1/4) * P_L + (1/4) * P_QR")
print("P_TC = (1/2) * P_SP + (1/2) * 0")
print("P_QR = (1/3) * P_L + (1/3) * P_KH + (1/3) * 1")
print("\nThis system can be written in matrix form Ax = b and solved for the probabilities.")

# The rooms are indexed as follows:
# 0: Main Entrance (ME), 1: Great Hall (GH), 2: Castle Kitchen (CK),
# 3: Secret Passage (SP), 4: Library (L), 5: Knights' Hall (KH),
# 6: Torture Chamber (TC), 7: Queen's Room (QR)

# We rearrange the equations into the standard form Ax = b.
# For example, P_ME - (1/2)P_GH - (1/2)P_CK = 0
A = np.array([
    [1.0, -1/2, -1/2, 0.0, 0.0, 0.0, 0.0, 0.0],
    [-1/4, 1.0, 0.0, -1/4, -1/4, -1/4, 0.0, 0.0],
    [-1/2, 0.0, 1.0, 0.0, 0.0, -1/2, 0.0, 0.0],
    [0.0, -1/4, 0.0, 1.0, 0.0, 0.0, -1/4, 0.0],
    [0.0, -1/3, 0.0, 0.0, 1.0, -1/3, 0.0, -1/3],
    [0.0, -1/4, -1/4, 0.0, -1/4, 1.0, 0.0, -1/4],
    [0.0, 0.0, 0.0, -1/2, 0.0, 0.0, 1.0, 0.0],
    [0.0, 0.0, 0.0, 0.0, -1/3, -1/3, 0.0, 1.0]
])

# The vector b contains the constant terms from the equations, which arise
# from direct paths to the Treasure Room or Vampire's Lair.
b = np.array([
    0.0,      # ME
    0.0,      # GH
    0.0,      # CK
    1/4,      # SP -> Treasure Room
    0.0,      # L
    0.0,      # KH
    0.0,      # TC
    1/3       # QR -> Treasure Room
])

# Solve the system of linear equations
probabilities = np.linalg.solve(A, b)

# The probability of success from the Main Entrance is the first element
prob_from_main_entrance = probabilities[0]

# Convert the floating-point result to an irreducible fraction
final_fraction = Fraction(prob_from_main_entrance).limit_denominator()

print("\nThe calculated probability of reaching the Treasure Room from the Main Entrance is:")
# The final equation is P(Success) = numerator / denominator
print(f"P(Treasure Room) = {final_fraction.numerator} / {final_fraction.denominator}")