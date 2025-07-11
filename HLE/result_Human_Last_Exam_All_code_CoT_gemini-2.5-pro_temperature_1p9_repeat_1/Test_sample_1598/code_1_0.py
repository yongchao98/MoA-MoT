# The number of leaf vertices in the K_1,n graph.
n = 100

# As determined by the analysis, the global labeling number for K_{1,n}
# is given by the formula 2^(n-1).
# We define the base and exponent for our calculation.
base = 2
exponent = n - 1

# Calculate the final result. Python handles large integers automatically.
result = base ** exponent

# The problem asks to output the numbers in the final equation.
# The final equation is: 2^99 = result.
print("The global labeling number of K_1,100 is given by the equation: base ^ (n - 1)")
print(f"n = {n}")
print(f"Base = {base}")
print(f"Exponent = {n} - 1 = {exponent}")
print(f"Result: {result}")
