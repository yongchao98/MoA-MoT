import math

# Calculate the value of alpha for the k=13 family
alpha_u = (13 + math.sqrt(165)) / 2

# Calculate the value of alpha for the k=5 family
alpha_v = (5 + math.sqrt(21)) / 2

# Calculate the limit L
L = 2 / math.log(alpha_u) + 2 / math.log(alpha_v)

# The question asks for the integer part of 10^4 * L
result = int(10000 * L)

print(result)