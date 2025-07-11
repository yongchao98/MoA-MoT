# Step 1: Define all X variables based on the puzzle's clues.
# From Clue 1 (Gesner, Kerosene, 1846): Y1 = 1846
X1, X2, X3, X4 = 1, 8, 4, 6

# From Clue 3 and structure of Y3 (Palladium-106 catalyst): Y3 = 4660
# Y3 = X3X4X8X6 => 46X86. So Y3 = 4660 gives X8=6, X6=0
X8, X6 = 6, 0

# From the hypothesis of the complete digit set {0..9} + one repeat,
# the remaining digits {2,3,5,7,9} are assigned to the remaining variables.
# Y4 = X9X10X11 is assigned {2,3,5} => Y4=235
X9, X10, X11 = 2, 3, 5
# Y2 = X5X68X7X6 is assigned {7,9} => Y2=70890
X5, X7 = 7, 9

# Step 2: Construct the full Y numbers.
Y1 = int(f"{X1}{X2}{X3}{X4}")
Y2 = int(f"{X5}{X6}{X2}{X7}{X6}")
Y3 = int(f"{X3}{X4}{X8}{X6}")
Y4 = int(f"{X9}{X10}{X11}")

# Step 3: Perform the final calculation, Y4 to the power of Y1.
result = Y4 ** Y1

# Step 4: Print the equation with all the numbers, as requested.
# The result is a very large number.
print(f"Based on the decoded clues:")
print(f"Y1 = {Y1}")
print(f"Y2 = {Y2}")
print(f"Y3 = {Y3}")
print(f"Y4 = {Y4}")
print("\nThe final calculation is Y4 to the power of Y1:")
print(f"{Y4} ** {Y1} = {result}")
