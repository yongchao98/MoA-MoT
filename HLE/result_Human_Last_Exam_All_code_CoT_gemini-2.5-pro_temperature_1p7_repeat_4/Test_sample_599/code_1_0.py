# Step 1: Define the position of the element we want to find.
n = 50

# Step 2: As determined by the analysis, the nth segmented number is given by the formula 2^(n-1).
# We can represent this as an equation: result = base ^ exponent.
base = 2
exponent = n - 1

# Step 3: Calculate the result.
result = base**exponent

# Step 4: Print the numbers in the final equation and the result.
# The final equation is 2^49 = result.
print(f"{base}^{exponent} = {result}")