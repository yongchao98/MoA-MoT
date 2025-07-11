X1 = 3
X2 = 2
X3 = 5
X4 = 10
X5 = 11
X6 = 9
X7 = 5
X8 = 1
X9 = 7

# The problem asks to print the final equation with each number.
print(f"{X1}'s {X2} + ({X3} X4 - {X5}'s {X6}) + {X7}'s {X4} = {X8}'s {X9}")
print(f"{X1} * {X2} + ({X3} * {X4} - {X5} * {X6}) + {X7} * {X4} = {X8} * {X9}")
result = X8 * X9
print(f"{X1*X2} + ({X3*X4} - {X5*X6}) + {X7*X4} = {result}")
print(f"{X1*X2} + {X3*X4 - X5*X6} + {X7*X4} = {result}")
print(f"{X1*X2 + (X3*X4 - X5*X6) + X7*X4} = {result}")
print("\nThe value of X12 is the result of this equation.")
print(f"X12 = {result}")
