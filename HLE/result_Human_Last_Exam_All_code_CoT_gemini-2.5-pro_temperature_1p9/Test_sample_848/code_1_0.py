import math

# The final result is derived from the formula:
# floor(10000 * (2 / ln((13 + sqrt(165)) / 2)))
# This python code calculates the value based on this formula.

# The numbers that appear in the final formula are:
num1 = 10000
num2 = 2
num3 = 13
num4 = 165
num5 = 2

# Per the instructions, we output each number in the final equation.
print(num1)
print(num2)
print(num3)
print(num4)
print(num5)

# Perform the calculation
phi = (num3 + math.sqrt(num4)) / num5
limit_L = num2 / math.log(phi)
result = num1 * limit_L

# The final answer is the integer part of the result.
final_answer = int(result)

# Output the final answer
print(final_answer)