# The numbers and operators are hidden in the colorful static noise in the image.
# By examining the area the character is pointing to, we can see the equation: 9 * 6 + 8 / 2

# Define the numbers found in the image
num1 = 9
num2 = 6
num3 = 8
num4 = 2

# Calculate the result following the order of operations (PEMDAS/BODMAS)
# Multiplication and division are performed before addition.
# (9 * 6) + (8 / 2) = 54 + 4 = 58
result = num1 * num2 + num3 / num4

# The prompt requires printing the equation showing each number.
# The result should be an integer.
print(f"Equation: {num1} * {num2} + {num3} / {num4} = {int(result)}")