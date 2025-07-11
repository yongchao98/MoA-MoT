import math

# Step 1: Set up the equation for the base of the asymptotic slice rank.
# Our calculated base of the slice rank growth. Let's call it `base_calc`.
base_calc = 2

# The given base of the slice rank growth from the problem. Let's call it `base_given`.
# The formula is (3 / (2**K))^n, so the base is 3 / (2**K).
# Let the numerator be `b` and the denominator part be `a**K`.
b = 3
a = 2

# Step 2: Formulate the equation to solve for K.
# The equation is base_calc = b / (a**K).
# In our case, this is 2 = 3 / (2**K).
print(f"The equation to solve is: {base_calc} = {b} / ({a}**K)")

# Step 3: Solve the equation for K.
# 2 * (2**K) = 3
# 2**(K+1) = 3
# K + 1 = log2(3)
# K = log2(3) - 1
# This is equivalent to K = log2(3/2).

K = math.log2(b / a)

# Step 4: Print the final value of K.
print(f"The value of K is: {K}")