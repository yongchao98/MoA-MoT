import math

print("To determine the most realistic upper bound O for a connected ε-envy-free allocation for 4 agents, we refer to recent findings in the field.")
print("A 2022 paper by Brânzei and Nisan provides a constructive algorithm with a query complexity of O(n^(n-1)).")
print("We will now calculate this bound for n=4 agents.")
print("-" * 25)

# Number of agents
n = 4

# The base of the exponentiation is n
base = n

# The exponent is n - 1
exponent = n - 1

# Calculate the result
result = int(math.pow(base, exponent))

print("The calculation is based on the formula: O = n ^ (n - 1)")
print(f"Substituting n = {n}:")
print(f"O = {base} ^ ({n} - 1)")
print(f"O = {base} ^ {exponent}")
print(f"O = {result}")
print("\nTherefore, the most realistic upper bound O that can be achieved is 64.")