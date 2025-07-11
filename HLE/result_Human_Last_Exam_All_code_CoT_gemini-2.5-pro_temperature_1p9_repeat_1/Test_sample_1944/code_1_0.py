# The obfuscated Javascript code uses specific patterns to create numbers.
# For example:
# 1 is created with `+!![]`
# 2 is created with `!![]+!![]`
# 3 is created with `!![]+!![]+!![]`

# Following the hint to form an equation, we can use these fundamental numbers.
num1 = 1
num2 = 2
num3 = 3

# We form a simple equation with these numbers.
result = num1 + num2 + num3

# As requested, we print the full equation, showing each number.
print(f"{num1} + {num2} + {num3} = {result}") 