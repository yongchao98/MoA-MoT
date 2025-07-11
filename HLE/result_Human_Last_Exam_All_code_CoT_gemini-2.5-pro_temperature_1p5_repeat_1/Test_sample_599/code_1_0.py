# The problem asks for the 50th segmented number.
# As determined by the analysis, the n-th segmented number is 2^(n-1).
n = 50

# The equation to find the 50th term is based on a power of 2.
base = 2
exponent = n - 1

# Calculate the result. Python handles large integers automatically.
result = base ** exponent

# The final equation is 2^49 = result.
# The instruction is to "output each number in the final equation".
# We will print the equation itself, which contains all the required numbers.
print(f"{base} ** {exponent} = {result}")
