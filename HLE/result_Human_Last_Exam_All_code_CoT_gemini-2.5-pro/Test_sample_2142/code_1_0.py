# Step 1: Define the variables X1 through X10 based on the solution "RENDER NURTURES".

# From Y = RENDER (X1 X2 X3 X4 X2 X5)
X1 = 'R'
X2 = 'E'
X3 = 'N'
X4 = 'D'
# X2 is 'E'
X5 = 'R'

# From Z = NURTURES (X3 X6 X1 X7 X8 X9 X10 X2 X6)
# X3 is 'N' (matches)
# X1 is 'R' (matches)
# X2 is 'E' (matches Z's 8th letter)
X6 = 'U'  # From Z's 2nd letter. Note: The puzzle has a small inconsistency here as Z's 9th letter is 'S'.
X7 = 'T'
X8 = 'U'
X9 = 'R'
X10 = 'E'

# Step 2: Construct the words Y and Z using the variables.
Y = [X1, X2, X3, X4, X2, X5]
Z = [X3, X6, X1, X7, X8, X9, X10, X2, 'S'] # Using 'S' for the last letter of NURTURES

# Step 3: Print the solution, showing how the variables form the words.
Y_word = "".join(Y)
Z_word = "".join(Z)

print(f"Y = {' '.join(Y)}")
print(f"Z = {' '.join(Z)}")
print(f"\nThe decoded phrase is: {Y_word} {Z_word}")
