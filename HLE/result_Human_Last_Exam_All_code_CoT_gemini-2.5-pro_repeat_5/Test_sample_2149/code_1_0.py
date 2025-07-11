import math

# Step 1: Assign the value to Y1 based on the analysis of Clue 1.
# The clue points to the year of the Hall-Héroult process discovery.
Y1 = 1886

# Step 2: Assign the value to Y4 based on the analysis of Clue 2.
# The clue points to a 3-digit number from the sequence in the constant 'e'.
Y4 = 459

# Step 3: Perform the calculation Y4 raised to the power of Y1.
# Python's ** operator can handle very large integers.
result = Y4 ** Y1

# Step 4: Print the equation and the final result.
# The problem requests to output each number in the final equation.
print(f"{Y4} ^ {Y1} = {result}")