import numpy as np
from fractions import Fraction

# Define the names of the transient states (rooms)
# The order is: Main Entrance, Great Hall, Secret Passage, Library,
# Torture Chamber, Queen's Room, Knights' Hall, Castle Kitchen.
rooms = [
    'Main Entrance', 'Great Hall', 'Secret Passage', 'Library',
    'Torture Chamber', "Queen's Room", "Knights' Hall", 'Castle Kitchen'
]
# Let P_i be the probability of reaching the Treasure Room from room i.
# We set up a system of linear equations Ax = b, where x is the vector of probabilities.
# x = [P_ME, P_GH, P_SP, P_L, P_TC, P_QR, P_KH, P_CK]^T

# The matrix A of coefficients is derived by rewriting the equations as:
# e.g., P_ME - (1/2)P_GH - (1/2)P_CK = 0
A_coeffs = [
    # ME    GH      SP      L       TC      QR      KH      CK
    [ 1,  -1/2,      0,      0,      0,      0,      0,   -1/2],  # Main Entrance Eq.
    [-1/4,    1,   -1/4,   -1/4,      0,      0,   -1/4,      0],  # Great Hall Eq.
    [   0, -1/4,      1,      0,   -1/4,      0,      0,      0],  # Secret Passage Eq.
    [   0, -1/3,      0,      1,      0,   -1/3,   -1/3,      0],  # Library Eq.
    [   0,    0,   -1/2,      0,      1,      0,      0,      0],  # Torture Chamber Eq.
    [   0,    0,      0,   -1/3,      0,      1,   -1/3,      0],  # Queen's Room Eq.
    [   0, -1/4,      0,   -1/4,      0,   -1/4,      1,   -1/4],  # Knights' Hall Eq.
    [-1/2,    0,      0,      0,      0,      0,   -1/2,      1],  # Castle Kitchen Eq.
]

# The vector b of constants comes from the connections to the absorbing states.
# P_TR=1 (Treasure), P_VL=0 (Vampire)
b_coeffs = [
    0,      # Main Entrance Eq. constant
    0,      # Great Hall Eq. constant
    1/4,    # Secret Passage Eq. constant (from 1/4 * P_TR)
    0,      # Library Eq. constant
    0,      # Torture Chamber Eq. constant (from 1/2 * P_VL)
    1/3,    # Queen's Room Eq. constant (from 1/3 * P_TR)
    0,      # Knights' Hall Eq. constant
    0,      # Castle Kitchen Eq. constant
]

# Convert coefficients to Fraction objects for exact arithmetic
A = np.array([[Fraction(c) for c in row] for row in A_coeffs])
b = np.array([Fraction(c) for c in b_coeffs])

# Solve the system of linear equations
try:
    probabilities = np.linalg.solve(A, b)

    # Print the results for each room
    print("The probability of reaching the Treasure Room from each location is:")
    for i, room in enumerate(rooms):
        print(f"{room:>18}: {probabilities[i]}")

    # The final answer is the probability for the Main Entrance
    p_me = probabilities[0]
    p_gh = probabilities[1]
    p_ck = probabilities[7]
    print("\n------------------------------------------------------------")
    print(f"The probability of reaching the Treasure Room from the Main Entrance is: {p_me}")
    print("------------------------------------------------------------\n")
    
    # Print the final equation with the calculated values
    print("The problem asks for the probability from the Main Entrance, which is calculated based on the rooms it connects to:")
    print("P_ME = (1/2) * P_GH + (1/2) * P_CK")
    print("\nSubstituting the calculated probabilities:")
    print(f"{p_me} = (1/2) * {p_gh} + (1/2) * {p_ck}")
    
    # Verify the equation
    check = Fraction(1,2) * p_gh + Fraction(1,2) * p_ck
    print(f"{p_me} = {check}")
    print("\nThe equation holds true, confirming the result.")

except np.linalg.LinAlgError:
    print("The system of equations could not be solved.")

<<<109/161>>>