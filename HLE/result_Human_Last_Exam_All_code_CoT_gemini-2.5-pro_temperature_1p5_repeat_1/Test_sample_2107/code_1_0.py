# Step 1: Define the variables X1 through X9 based on the riddle's clues.
X1 = 8
X2 = 2
X3 = 5
X4 = 7
X5 = 10
X6 = 5
X7 = 14
X8 = 11
X9 = 9

# Step 2: Verify the first mathematical equation.
# The equation is: X1 * X2 + (X3 * X4 - X5 * X6) + X7 * X4 = X8 * X9
left_hand_side = X1 * X2 + (X3 * X4 - X5 * X6) + X7 * X4
right_hand_side = X8 * X9

# Step 3: Print the equation with the solved numbers.
print("Verifying the equation:")
print(f"{X1} * {X2} + ({X3} * {X4} - {X5} * {X6}) + {X7} * {X4} = {X8} * {X9}")
print(f"{left_hand_side} = {right_hand_side}")

if left_hand_side == right_hand_side:
    print("The equation holds true.")
else:
    print("There is a mismatch in the equation.")

# Step 4: Determine the final unknown, X12.
# Based on the analysis of the second line: X8's X9's X10 X11 X12.
# X10 is derived from the Prutkov clue, where the final word "братьев" (brothers) has 7 letters.
# The puzzle's structure suggests X12 is revealed by a previous clue in the sequence.
# Thus, X12 takes the value of X10.
X10 = 7
X12 = 7

# Step 5: Print the final answer for X12.
print("\nThe solution for X12 is:")
print(f"X12 = {X12}")

<<<7>>>