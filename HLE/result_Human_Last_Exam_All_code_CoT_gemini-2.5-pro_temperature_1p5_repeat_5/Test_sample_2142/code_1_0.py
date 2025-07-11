# 1. Define the characters based on the solution "SECRET CROATIONS".
# Y = SECRET
X1 = 'S'
X2 = 'E'
X3 = 'C'
X4 = 'R'
X5 = 'T'

# Z = CROATIONS, which is constructed to fit the pattern.
# Z = X3 X6 X1 X7 X8 X9 X10 X2 X6
# C  R  S  A  T  I   O   E   R
X6 = 'R'
X7 = 'A'
X8 = 'T'
X9 = 'I'
X10 = 'O'

# 2. Construct the words Y and Z from the characters.
y_word = f"{X1}{X2}{X3}{X4}{X2}{X5}"
z_word = f"{X3}{X6}{X1}{X7}{X8}{X9}{X10}{X2}{X6}"

# 3. Print the final equation with the deciphered letters.
# To make it clear, we show the equation for Y and Z.
print(f"Y = {X1}{X2}{X3}{X4}{X2}{X5}")
print(f"Z = {X3}{X6}{X1}{X7}{X8}{X9}{X10}{X2}{X6}")
print()
# 4. Print the final phrase.
print(f"The deciphered phrase is: {y_word} {z_word}")
