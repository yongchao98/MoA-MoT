# Step 1: Define the variables X1 through X10 based on the decoded phrase "CUTOUT TECHNIQUE".
X1 = 'C'
X2 = 'U'
X3 = 'T'
X4 = 'O'
# X5 is also 'T', the same as X3.
X5 = 'T' 
X6 = 'E'
X7 = 'H'
X8 = 'N'
X9 = 'I'
X10 = 'Q'

# Step 2: Construct the words Y and Z using the given formulas.
# Y = X1X2X3X4X2X5
Y = f"{X1}{X2}{X3}{X4}{X2}{X5}"

# Z = X3X6X1X7X8X9X10X2X6
Z = f"{X3}{X6}{X1}{X7}{X8}{X9}{X10}{X2}{X6}"

# Step 3: Print the final deciphered phrase "Y Z".
# The prompt asks to output each component of the final equation, which are the words Y and Z.
print(f"{Y} {Z}")