import numpy as np

# Step 1: Define the constants lambda1 and lambda2
sqrt17 = np.sqrt(17)
lambda1 = (1 + sqrt17) / 2
lambda2 = (1 - sqrt17) / 2

# Step 2: State the assumption for the unknown initial condition
x2_0 = 0

# Step 3: Calculate the numerical components of the expression
# The expression is: P1 * x2(0) - P2 - P3
# Component 1: The coefficient of x2(0)
P1 = (2/3 * lambda1 * np.exp(lambda1 / 2)) - (1/3 * np.exp(lambda1 / 2))

# Component 2: The second term to be subtracted
P2 = (2/3 * lambda2 * np.exp(lambda2 / 2))

# Component 3: The third term to be subtracted
P3 = (10/3 * np.exp(lambda1 / 2))

# Step 4: Calculate the final result using the assumption
result = P1 * x2_0 - P2 - P3

# Step 5: Print the components of the final equation and the result
print("The expression to evaluate is:")
print(f"({P1}) * x2(0) - ({P2}) - ({P3})")
print(f"\nAssuming the initial condition x2(0) = {x2_0}, the expression becomes:")
print(f"({P1}) * {x2_0} - ({P2}) - ({P3})")
print(f"= {P1 * x2_0} - ({P2}) - ({P3})")
print(f"= {result}")

# Final numerical result
print(f"\nThe final calculated value is: {result}")