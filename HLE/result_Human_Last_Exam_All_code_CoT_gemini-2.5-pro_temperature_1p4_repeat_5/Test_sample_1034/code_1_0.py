import math

# Let S = SNYWR. We will assume a 5-digit integer value for S.
S = 98765

# In triangle ABC, we are given:
# Side BC = S
# Angle B = 30 degrees
# We are looking for the number of integer values of side AC for which
# there are two possible values for side AB.

# This is the ambiguous case of the Law of Sines. For two distinct triangles
# to be possible, the length of side AC must be greater than the altitude
# from vertex C but less than the length of side BC.

# Altitude h = BC * sin(B) = S * sin(30) = S / 2
# The condition is: S / 2 < AC < S

# We need to count the number of integers for AC in this range.
lower_bound = S / 2
upper_bound = S

# The smallest integer for AC is floor(S / 2) + 1
first_integer_ac = math.floor(lower_bound) + 1

# The largest integer for AC is S - 1
last_integer_ac = upper_bound - 1

# The number of integer values is the count of integers in the range
# [first_integer_ac, last_integer_ac], inclusive.
# The final equation to find the count is: (last - first + 1)
count = last_integer_ac - first_integer_ac + 1

print(f"Given S = {S}, the condition for AC is:")
print(f"{S} / 2 < AC < {S}")
print(f"Which means: {lower_bound} < AC < {upper_bound}")
print("\nWe need to count the integers in this range.")
print(f"The first valid integer for AC is {first_integer_ac}.")
print(f"The last valid integer for AC is {last_integer_ac}.")
print("\nThe final equation for the count is: Last Integer - First Integer + 1")
print(f"So, the number of integer values is: {last_integer_ac} - {first_integer_ac} + 1 = {count}")
