# Step 1: Define the values for each character based on the decoded words
# Y = CUTUPS, Z = TECHNIQUE
X1 = 'C'
X2 = 'U'
X3 = 'T'
X4 = 'P'
X5 = 'S'
X6 = 'E'
X7 = 'H'
X8 = 'N'
X9 = 'I'
X10 = 'Q'

# Step 2: Create lists of characters for Y and Z based on the problem's equations
Y_chars = [X1, X2, X3, X4, X2, X5]
Z_chars = [X3, X6, X1, X7, X8, X9, X10, X2, X6]

# Step 3: Print the characters of each word in the final equation, separated by spaces.
# The '*' unpacks the lists, and the print function places a space between each argument.
# An extra ' ' is added to create a wider space between the two words.
print(*Y_chars, ' ', *Z_chars)
