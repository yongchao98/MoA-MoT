# The Piet image's structure suggests a mathematical equation.
# The image is divided into two sections by a black bar, representing two numbers.

# In the first section, a 6-color cycle repeats 2 times.
# This gives us the first number.
num1 = 2

# In the second section, another 6-color cycle also repeats 2 times.
# This gives us the second number.
num2 = 2

# We assume the operation is addition, the most fundamental one.
result = num1 + num2

# The program should print the full equation, including each number and the result.
print(f"{num1} + {num2} = {result}")