# The first number is derived from the first color on the border, Red.
# In the chosen color-to-number mapping (based on RGB bit values), Red corresponds to 4.
n1 = 4

# The second number is from the second color, Green, which maps to 2.
n2 = 2

# The operation is assumed to be addition, suggested by the additive color model (Red + Green = Yellow).
# The result of the addition is 4 + 2 = 6.
result = n1 + n2

# The number 6 maps to the color Yellow, which is present in the border sequence, confirming the equation.
# The code now prints the final equation, showing each number.
print(f"{n1} + {n2} = {result}")