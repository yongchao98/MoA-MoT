# Step 1 & 2: Decode Clue 1 to find Y1 and the first set of X variables.
# The clue refers to the Hall-Héroult process, discovered in 1886.
X1 = 1
X2 = 8
X3 = 8
X4 = 6

# Construct Y1 from its component digits.
Y1 = int(f"{X1}{X2}{X3}{X4}")

# Step 3 & 4: Decode Clue 2 to find the remaining X variables.
# The clue refers to a sequence of numbers (24, 735, 987) from
# Georges Perec's "Life: A User's Manual". We take the unique digits in order.
X5 = 2
X6 = 4
X7 = 7
X8 = 3
X9 = 5
X10 = 9
X11 = 8

# Step 6: Calculate Y4 using its component X variables.
Y4 = int(f"{X9}{X10}{X11}")

# Step 7: Formulate and print the final answer as an equation.
# The prompt asks for "Y4 to the Y1", which suggests an exponentiation.
# As instructed, we output each number in the final equation.

print("The value for Y4 has been calculated as:")
print(Y4)
print("The value for Y1 has been calculated as:")
print(Y1)
print("\nThe final equation for the Y4 to the Y1-Hall topological state indices is:")
# This prints the final equation with the computed numbers.
print(f"{Y4} ** {Y1}")
