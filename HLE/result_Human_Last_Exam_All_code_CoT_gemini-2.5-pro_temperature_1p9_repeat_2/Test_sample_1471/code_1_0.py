# 1. Define the values based on interpretation of the numerals.
# ⅣⅭⅯ -> 4 + 100 + 1000
num1 = 1104
# 𒌋𒌋𒌋 -> 10 + 10 + 10
num2 = 30
# 𒌋𒌋𒁹𒁹𒁹 -> 10 + 10 + 1 + 1 + 1
num3 = 23

# 2. Perform the calculation following order of operations (division first).
result = (num1 / num2) - num3

# 3. Print the equation with the final answer.
# Using an f-string to format the output nicely.
print(f"{num1} / {num2} - {num3} = {result}")