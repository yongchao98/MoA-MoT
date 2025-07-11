# Step 1: Define the characters X1 through X10 based on the words GENDER and NEGATIONS.
# From Y = GENDER
X1 = 'G'
X2 = 'E'
X3 = 'N'
X4 = 'D'
X5 = 'R'

# From Z = NEGATIONS
# Note: Z's formula (X3X6X1...X2X6) requires Z[1]==Z[8]. In 'NEGATIONS', Z[1]=='E' and Z[8]=='S', which doesn't match.
# We proceed by assigning the letters from NEGATIONS to the remaining X variables to fulfill the task.
X6 = 'E'  # Second letter of NEGATIONS
X7 = 'A'  # Fourth letter
X8 = 'T'  # Fifth letter
X9 = 'I'  # Sixth letter
X10 = 'O' # Seventh letter

# Step 2: Build the words Y and Z using the formulas from the riddle.
Y = X1 + X2 + X3 + X4 + X2 + X5
Z_formula = [X3, X6, X1, X7, X8, X9, X10, X2, X6]
Z = "".join(Z_formula)

# Step 3: Print the final deciphered phrase "Y Z" and show the work.
print(f"Y = {X1}+{X2}+{X3}+{X4}+{X2}+{X5} = {Y}")
print(f"Z = {X3}+{X6}+{X1}+{X7}+{X8}+{X9}+{X10}+{X2}+{X6} = {Z}")
print("")
print(f"The deciphered phrase is: {Y} {Z}")
