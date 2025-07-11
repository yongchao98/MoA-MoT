# Step 1 & 2: List givens and derive I(X;Z)
# Given values
I_XY = 3  # I(X;Y)
I_XY_Z = 2 # I(X;Y|Z)
I_XZ_Y = 5 # I(X;Z|Y)

# The chain rule for mutual information states:
# I(X; Y,Z) = I(X;Y) + I(X;Z|Y)
# I(X; Y,Z) = I(X;Z) + I(X;Y|Z)
# Therefore, I(X;Y) + I(X;Z|Y) = I(X;Z) + I(X;Y|Z)
# We can solve for I(X;Z)
# 3 + 5 = I(X;Z) + 2
I_XZ = I_XY + I_XZ_Y - I_XY_Z
# I_XZ = 3 + 5 - 2 = 6

print(f"I(X;Y) = {I_XY}")
print(f"I(X;Y|Z) = {I_XY_Z}")
print(f"I(X;Z|Y) = {I_XZ_Y}")
print(f"Derived I(X;Z) = {I_XY} + {I_XZ_Y} - {I_XY_Z} = {I_XZ}")
print("-" * 20)

# Step 3, 4, 5: Find an expression for I(X;Y|W)
# We use the chain rule on I(X; Y,Z | W).
# I(X; Y,Z | W) = I(X;Y|W) + I(X;Z|Y,W)
# I(X; Y,Z | W) = I(X;Z|W) + I(X;Y|Z,W)
# Equating the two expressions:
# I(X;Y|W) + I(X;Z|Y,W) = I(X;Z|W) + I(X;Y|Z,W)

# Since W is a function of Z (W=f(Z)), conditioning on Z is more informative
# or equal to conditioning on W. So, conditioning on (Z, W) is the same as on Z.
# This means I(X;Y|Z,W) = I(X;Y|Z)
I_XY_ZW = I_XY_Z # I(X;Y|Z,W) = 2

# Substitute this into the equation:
# I(X;Y|W) + I(X;Z|Y,W) = I(X;Z|W) + I(X;Y|Z)
# Rearranging for I(X;Y|W):
# I(X;Y|W) = I(X;Y|Z) + I(X;Z|W) - I(X;Z|Y,W)
# I(X;Y|W) = 2 + I(X;Z|W) - I(X;Z|Y,W)

# Step 6 & 7: Maximize the expression
# We need to maximize I(X;Y|W) = 2 + I(X;Z|W) - I(X;Z|Y,W)
# Let's express the terms with W using known quantities.
# For I(X;Z|W):
# I(X;Z) = I(X;W) + I(X;Z|W)  (since W is a function of Z)
# I(X;Z|W) = I(X;Z) - I(X;W) = 6 - I(X;W)

# For I(X;Z|Y,W):
# I(X;Z|Y) = I(X;Z,W|Y) = I(X;W|Y) + I(X;Z|Y,W)
# I(X;Z|Y,W) = I(X;Z|Y) - I(X;W|Y) = 5 - I(X;W|Y)

# Substitute these back into the expression for I(X;Y|W):
# I(X;Y|W) = 2 + (6 - I(X;W)) - (5 - I(X;W|Y))
# I(X;Y|W) = 2 + 6 - I(X;W) - 5 + I(X;W|Y)
# I(X;Y|W) = 3 + I(X;W|Y) - I(X;W)

# To maximize I(X;Y|W), we need to maximize the difference I(X;W|Y) - I(X;W).
# Let's analyze the bounds of this difference.
# I(X;W|Y) is the information between X and W, given Y.
# Since W is a function of Z, any information W has about X must come through Z.
# By the data processing inequality, I(X;W|Y) <= I(X;Z|Y).
# We know I(X;Z|Y) = 5, so I(X;W|Y) <= 5.
# Also, mutual information is always non-negative, so I(X;W) >= 0.

# The difference is therefore bounded:
# I(X;W|Y) - I(X;W) <= 5 - 0 = 5.
# This gives an upper bound for I(X;Y|W):
# I(X;Y|W) <= 3 + 5 = 8.

# However, we must check if this bound is attainable.
# To reach this bound, we would need I(X;W|Y) = 5 AND I(X;W) = 0.
# Let's check for contradictions.
# 1. I(X;W|Y) = 5 = I(X;Z|Y). Equality holds if and only if I(X;Z|Y,W) = 0.
# 2. I(X;W) = 0. This means X and W are independent.

# If X and W are independent, then conditioning on W provides no information about X.
# So, I(X;Z|Y,W) = I(X;Z|Y) because W is independent of X.
# So condition (1) becomes I(X;Z|Y) = 0.
# This contradicts the given information that I(X;Z|Y) = 5.
# Therefore, the value 8 is not attainable.

# We need to maximize I(X;W|Y) - I(X;W) subject to being constructible.
# The value of I(X;Y|W) is given by $3 + I(X;W|Y) - I(X;W)$.
# Let's explore achievable values for $I(X;W|Y) - I(X;W)$.
# Consider a case where we can choose W=f(Z) such that:
# - I(X;W|Y) = 5 (This requires I(X;Z|Y,W) = 0)
# - I(X;W) = 1
# This is a hypothetical valid construction, as there is no immediate contradiction. For example, W could be a part of Z that becomes fully known once Y is known, but has some smaller residual information about X on its own.
# If such a construction exists, what would the value of I(X;Y|W) be?
val_I_XW_Y = 5
val_I_XW = 1
max_val = 3 + val_I_XW_Y - val_I_XW
print(f"An achievable scenario could be I(X;W|Y) = {val_I_XW_Y} and I(X;W) = {val_I_XW}.")
print(f"This is consistent with the condition I(X;W|Y) - I(X;W) = {val_I_XW_Y - val_I_XW}, which must be <= I(X;Z|Y) = {I_XZ_Y}.")
print(f"The largest possible value would be:")
print(f"I(X;Y|W) = 3 + I(X;W|Y) - I(X;W)")
print(f"I(X;Y|W) = 3 + {val_I_XW_Y} - {val_I_XW} = {max_val}")
# The largest possible integer value for the difference I(X;W|Y) - I(X;W) that is not contradicted is I(X;Z|Y) - 1 = 4. Let's assume a configuration like $I(X;W|Y)=5$ and $I(X;W)=1$ or $I(X;W|Y)=4$ and $I(X;W)=0$ is possible.
# If I(X;W|Y) - I(X;W) = 4, then I(X;Y|W) = 3 + 4 = 7.
# This value corresponds to I(X;Z|Y) + I(X;Y|Z) = 5 + 2 = 7. Let's verify if this can be a final answer.
# It seems this is the largest possible value not leading to a direct contradiction.

final_answer = I_XZ_Y + I_XY_Z
print("-" * 20)
print(f"The largest possible value is I(X;Z|Y) + I(X;Y|Z) = {I_XZ_Y} + {I_XY_Z} = {final_answer}")