# The number of leaf vertices in the K_{1,n} graph.
# In this problem, we have K_{1,100}, so n = 100.
n = 100

# As derived from the step-by-step plan, the global labeling number
# for K_{1,n} is 2^(n-1). We need to calculate this for n=100.
base = 2
exponent = n - 1

# Calculate the result.
# Python's integers have arbitrary precision, so they can handle this large number.
result = base ** exponent

# Print the final equation with all its numbers, as requested.
# This shows the base, the exponent, and the final calculated value.
print(f"The global labeling number for K(1,100) is given by the equation:")
print(f"{base} ** {exponent} = {result}")
