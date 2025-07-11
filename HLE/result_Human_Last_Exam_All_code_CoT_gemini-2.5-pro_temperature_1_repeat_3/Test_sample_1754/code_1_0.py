import math

# Step 1: Define the given constants lambda_1 and lambda_2
lambda1 = (1 + math.sqrt(17)) / 2
lambda2 = (1 - math.sqrt(17)) / 2

# Step 2: Assume x2(0) = 0 based on the context of the controllability problem starting from the origin.
x2_0 = 0

# Step 3: Calculate each numerical part of the equation as requested by the prompt.
# The equation is of the form: C1 * x2(0) - C2 - C3

# C1 = (2/3 * lambda1 * exp(lambda1/2) - 1/3 * exp(lambda1/2))
C1 = (2/3) * lambda1 * math.exp(lambda1 / 2) - (1/3) * math.exp(lambda1 / 2)

# C2 = 2/3 * lambda2 * exp(lambda2/2)
C2 = (2/3) * lambda2 * math.exp(lambda2 / 2)

# C3 = 10/3 * exp(lambda1/2)
C3 = (10/3) * math.exp(lambda1 / 2)

# Step 4: Calculate the final result
result = C1 * x2_0 - C2 - C3

# Step 5: Print the final equation with the computed values
print("The expression to calculate is: (C1) * x2(0) - (C2) - (C3)")
print(f"Assuming x2(0) = {x2_0}, the calculated values are:")
print(f"C1 = {C1}")
print(f"C2 = {C2}")
print(f"C3 = {C3}")
print("\nThe final equation is:")
print(f"({C1}) * {x2_0} - ({C2}) - ({C3}) = {result}")
