import numpy as np

# --- 1. Define the transition rates ---
l01 = 0.019
l10 = 0.65
l12 = 0.4
l21 = 0.392
l23 = 0.008
l31 = 0.008

# --- 2. Set up the system of linear equations A*p = b for the steady state ---
# P = [P0, P1, P2, P3]
# The equations are:
# l01*P0 - l10*P1 = 0
# l12*P1 - (l21+l23)*P2 = 0
# l23*P2 - l31*P3 = 0
# P0 + P1 + P2 + P3 = 1

A = np.array([
    [l01, -l10, 0, 0],
    [0, l12, -(l21 + l23), 0],
    [0, 0, l23, -l31],
    [1, 1, 1, 1]
])

b = np.array([0, 0, 0, 1])

# --- 3. Solve for the probabilities P ---
try:
    P = np.linalg.solve(A, b)
    P0_inf = P[0]
    P1_inf = P[1]

    # --- 4. Calculate the required sum ---
    result = P0_inf + P1_inf

    # --- 5. Print the results in a detailed format ---
    # Remember in the final code you still need to output each number in the final equation!
    print(f"The steady-state probabilities are:")
    print(f"P0(inf) = {P0_inf}")
    print(f"P1(inf) = {P1_inf}")
    print(f"P2(inf) = {P[2]}")
    print(f"P3(inf) = {P[3]}")
    print(f"\nThe required sum is P0(inf) + P1(inf):")
    print(f"{P0_inf} + {P1_inf} = {result}")

except np.linalg.LinAlgError:
    print("The system of equations is singular and cannot be solved.")
