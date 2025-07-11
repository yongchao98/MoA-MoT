# Step 1: Decode the numbers from the equation.
# The first number is represented by four Egyptian hieroglyphs for 100,000.
num1 = 400000

# The second number is in Babylonian cuneiform: 𒌋𒌋𒌋 (10 + 10 + 10).
num2 = 30

# The third number is also in Babylonian cuneiform: 𒌋𒌋𒁹𒁹𒁹 (10 + 10 + 1 + 1 + 1).
num3 = 23

# Step 2: Perform the calculation.
# Following the order of operations (division then subtraction).
result = num1 / num2 - num3

# Step 3: Print the equation with modern numbers and the result.
print(f"The equation with modern numbers is:")
print(f"{num1} / {num2} - {num3} = {result}")
