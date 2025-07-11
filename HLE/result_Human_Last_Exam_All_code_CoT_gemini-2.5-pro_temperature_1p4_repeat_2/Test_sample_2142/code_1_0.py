# Define the characters based on the deciphered words.
X1 = 'G'
X2 = 'E'
X3 = 'N'
X4 = 'D'
X5 = 'R'
X6 = 'E'
X7 = 'A'
X8 = 'T'
X9 = 'I'
X10 = 'V'

# Define the words Y and Z using the characters.
Y = X1 + X2 + X3 + X4 + X2 + X5
Z = X3 + X6 + X1 + X7 + X8 + X9 + X10 + X2 + X6

# Print the step-by-step solution as per the instructions.
# "output each number in the final equation!"
# I interpret this as showing the substitution for each position.

print("Decoding Y = X1X2X3X4X2X5:")
print(f"X1={X1}, X2={X2}, X3={X3}, X4={X4}, X2={X2}, X5={X5}")
print(f"Resulting Y: {Y}\n")

print("Decoding Z = X3X6X1X7X8X9X10X2X6:")
print(f"X3={X3}, X6={X6}, X1={X1}, X7={X7}, X8={X8}, X9={X9}, X10={X10}, X2={X2}, X6={X6}")
print(f"Resulting Z: {Z}\n")

print("The decoded phrase 'Y Z' is:")
print(f"{Y} {Z}")