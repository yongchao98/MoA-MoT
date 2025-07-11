# Step 1: Assign the determined character values to variables X1 through X10.
X1 = 'G'
X2 = 'E'
X3 = 'N'
X4 = 'D'
X5 = 'R'
X6 = 'E' # Solved from Z
X7 = 'O'
X8 = 'T'
X9 = 'I'
X10 = 'A'

# Step 2: Construct the words Y and Z using the variables as per the riddle's formulas.
Y = X1 + X2 + X3 + X4 + X2 + X5
Z = X3 + X6 + X1 + X7 + X8 + X9 + X10 + X2 + X6

# Step 3: Print the breakdown of each word and the final result.
# The request was to "output each number in the final equation". 
# This is interpreted as showing each character component of the final words.

print("The deciphered words are Y and Z.")
print("Solving for Y = X1X2X3X4X2X5:")
print(f"Y = {X1} {X2} {X3} {X4} {X2} {X5}")
print(f"Y = {Y}")

print("\nSolving for Z = X3X6X1X7X8X9X10X2X6:")
print(f"Z = {X3} {X6} {X1} {X7} {X8} {X9} {X10} {X2} {X6}")
print(f"Z = {Z}")

print("\nP-Orridge's deep commitment was to Y Z:")
print(f"{Y} {Z}")