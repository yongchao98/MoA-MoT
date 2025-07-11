# The number of leaf vertices in the K_1,n graph is n.
n = 100

# As derived in the explanation, the global labeling number for K_1,n is 2^(n-1).
# We need to calculate this for n=100.
base = 2
exponent = n - 1

# Calculate the result
result = base ** exponent

# The problem asks to output each number in the final equation.
# Here we print the base, exponent, and the final calculated value.
print(f"The global labeling number of K(1,100) is given by the equation: {base}^{exponent}")
print(f"The result is: {result}")
