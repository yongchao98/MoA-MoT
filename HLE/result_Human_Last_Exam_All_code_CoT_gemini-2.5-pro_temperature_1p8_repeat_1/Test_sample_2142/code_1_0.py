# Step 1: Define the mapping from the variables X1-X10 to the correct letters.
# This mapping was determined by solving the cryptogram based on the puzzle's constraints.
mapping = {
    'X1': 'E',
    'X2': 'N',
    'X3': 'T',
    'X4': 'I',
    'X5': 'Y',
    'X6': 'C',
    'X7': 'O',
    'X8': 'R',
    'X9': 'L',
    'X10': 'D'
}

# Step 2: Define the structure of Y and Z using the variable names.
y_formula = ['X1', 'X2', 'X3', 'X4', 'X2', 'X5']
z_formula = ['X3', 'X6', 'X1', 'X7', 'X8', 'X9', 'X10', 'X2', 'X6']

# Step 3: Decode the words Y and Z by substituting the variables with their mapped letters.
Y = "".join([mapping[var] for var in y_formula])
Z = "".join([mapping[var] for var in z_formula])

# Step 4: Print the final output, showing the decoded words.
# The prompt requires showing the decoded words that form the answer.
print(f"Y = {Y}")
print(f"Z = {Z}")
print(f"The deciphered phrase is: {Y} {Z}")
