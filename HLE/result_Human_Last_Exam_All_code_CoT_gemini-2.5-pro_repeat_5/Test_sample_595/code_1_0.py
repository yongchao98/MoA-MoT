import math

# Define properties of the triangle and its optimal placement
side_leg = 18
# In the optimal orientation, the legs are broken down into x and y components.
# The component size 'd' is such that d*sqrt(2) = 18.
d = 9 * math.sqrt(2)

# The number of grid lines crossed depends on the integer part of d.
I = math.floor(d)

# Explanation of the calculation for the maximum number of crossings (k).
# For a line segment, the number of squares crossed is the sum of horizontal and vertical grid lines crossed.
# By carefully choosing the fractional part of the triangle's starting vertex, we can maximize this for each side.

# For a segment with projection 'd' on an axis, the number of crossings can be floor(d) or floor(d)+1. The max is I+1.
# For a segment with projection '2d', the max number of crossings is 2*I+2.

# 1. First leg (length 18):
# Projections on x and y axes are both 'd'.
# Maximum crossings = (I + 1) + (I + 1)
k_leg1 = (I + 1) + (I + 1)

# 2. Second leg (length 18):
# Projections on x and y axes are also 'd'.
# Maximum crossings = (I + 1) + (I + 1)
k_leg2 = (I + 1) + (I + 1)

# 3. Hypotenuse (length 18*sqrt(2)):
# This side is aligned with an axis, so its projection on that axis is 2d = 18*sqrt(2).
# Maximum crossings = 2*I + 2
k_hypotenuse = 2 * I + 2

# The largest number k is the sum of the maximums for each side.
k_total = k_leg1 + k_leg2 + k_hypotenuse

print(f"The side lengths of the triangle are {side_leg}, {side_leg}, and {side_leg}*sqrt(2).")
print("To maximize the number of grid squares crossed, the triangle's legs are oriented at 45 degrees to the grid axes.")
print(f"The number of squares crossed by the first leg is {k_leg1}.")
print(f"The number of squares crossed by the second leg is {k_leg2}.")
print(f"The number of squares crossed by the hypotenuse is {k_hypotenuse}.")
print(f"The largest number k is the sum of these values:")
print(f"{k_leg1} + {k_leg2} + {k_hypotenuse} = {k_total}")
