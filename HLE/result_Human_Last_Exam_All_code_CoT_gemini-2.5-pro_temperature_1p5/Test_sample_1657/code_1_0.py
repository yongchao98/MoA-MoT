# Based on the analysis, the problem is unsolvable for two main reasons:
# 1. Missing Information: The initial distance to Pandora is a required input to
#    calculate the travel time from Earth's frame of reference (a), but it is not provided.
# 2. Computational Limitation: Calculating the onboard time (b) requires finding a
#    square root for time dilation. The Wuxing architecture's 'frac' data type has
#    a 'signed char' numerator, which is too small to handle the intermediate values
#    needed to implement a square root algorithm, making the calculation impossible.

# As per the problem's instructions, if a program cannot be written to solve for
# a and b, the answer should be 0:0.

a = 0
b = 0

# The problem asks for the answer in the format a:b.
print(f"{a}:{b}")