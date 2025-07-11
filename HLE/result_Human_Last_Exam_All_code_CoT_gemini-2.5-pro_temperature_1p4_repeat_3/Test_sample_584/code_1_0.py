# Parameters from the problem statement
r = 3
b = 9
false_positive_rate = 0.05

# The equation for the threshold 's' is derived from the S-curve formula:
# P(s) = 1 - (1 - s^r)^b
# We need to find 's' when P(s) = false_positive_rate.
#
# So, we solve for 's' in:
# false_positive_rate = 1 - (1 - s^r)^b
# s = (1 - (1 - false_positive_rate)^(1/b))^(1/r)

# Perform the calculation step-by-step
# Step 1: 1 - false_positive_rate
val1 = 1 - false_positive_rate

# Step 2: (val1)^(1/b)
val2 = val1**(1/b)

# Step 3: 1 - val2
val3 = 1 - val2

# Step 4: (val3)^(1/r)
s = val3**(1/r)

# Round the final result to three decimal points
s_rounded = round(s, 3)

# Print the calculation process showing each number in the equation
print("The calculation for the threshold 's' is as follows:")
print(f"s = (1 - (1 - {false_positive_rate})^(1/{b}))^(1/{r})")
print(f"s = (1 - ({val1})^(1/{b}))^(1/{r})")
print(f"s = (1 - {val2:.5f})^(1/{r})")
print(f"s = ({val3:.5f})^(1/{r})")
print(f"s = {s:.3f}")
print("\nWhen rounded to three decimal points, the threshold is:")
print(s_rounded)
