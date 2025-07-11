# Assign the decoded numerical values to variables.
X1 = 10
X2 = 2
X3 = 50
X4 = 42
X5 = 1000000
X6 = 25
X7 = 42
X8 = 11
X9 = 6

# The first equation is a red herring, but per the instructions,
# here are the numbers substituted into it.
# The operations ('s and space) are interpreted as multiplication for display.
equation_str = f"{X1} * {X2} + ({X3} * {X4} - {X5} * {X6}) + {X7} * {X4} = {X8} * {X9}"

print("The equation from the riddle with the decoded numbers is:")
print(equation_str)
print("\nThis equation is a red herring, as it does not hold true (-24996116 != 66).")

# The actual solution comes from the logical deduction that the theme is "first".
X12 = 1

print("\nThe value of X12, based on the logical puzzle in the second line, is:")
print(X12)