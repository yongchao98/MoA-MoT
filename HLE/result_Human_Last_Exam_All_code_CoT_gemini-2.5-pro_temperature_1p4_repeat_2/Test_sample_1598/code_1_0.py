# The graph is K_1,n where n=100.
n = 100

# As derived in the explanation, the global labeling number for K_1,n is 2^(n-1).
# The final equation is: k = base ^ exponent
# We will define these numbers and calculate the result.

base = 2
exponent = n - 1

# Calculate the final result using Python's support for large integers.
result = base ** exponent

# Print the components of the final equation and the result, as requested.
print("The global labeling number of the graph K_1,100 is determined by the equation:")
print(f"  k = {base}^{exponent}")
print("\nEach number in this final equation is:")
print(f"  base = {base}")
print(f"  exponent = {exponent}")
print(f"\nThe resulting global labeling number k is:")
print(result)