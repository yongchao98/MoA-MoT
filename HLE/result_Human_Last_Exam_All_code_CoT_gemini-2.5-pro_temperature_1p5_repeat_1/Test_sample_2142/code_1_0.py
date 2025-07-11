# 1. Decipher the characters X1 through X10 based on the thematic answer.
# The phrase is "ENGLISH HERETIC", which contains 10 unique letters.
X1 = 'E'
X2 = 'N'
X3 = 'G'
X4 = 'L'
X5 = 'I'
X6 = 'S'
X7 = 'H'
X8 = 'R'
X9 = 'T'
X10 = 'C'

# 2. Construct the words Y and Z using the formulas from the prompt.
Y_list = [X1, X2, X3, X4, X2, X5]
Z_list = [X3, X6, X1, X7, X8, X9, X10, X2, X6]

Y = "".join(Y_list)
Z = "".join(Z_list)

# 3. Print the final phrase and the equations that form it.
print(f"The deciphered phrase is: {Y} {Z}")
print("\nThis is constructed as follows:")
print(f"Y = {' '.join(Y_list)}")
print(f"Z = {' '.join(Z_list)}")
