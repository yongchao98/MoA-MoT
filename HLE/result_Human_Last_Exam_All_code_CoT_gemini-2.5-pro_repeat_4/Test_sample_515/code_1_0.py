import math

# The problem is to find the probability that a 2D simple random walk
# starting at (0,1) and conditioned to never enter the origin
# will eventually hit the set of neighbors of (3600, 0).

# This probability is well-approximated by the "escape probability" for a
# standard simple random walk: the probability that a walk starting at a neighbor
# of the origin reaches a circle of radius R before returning to the origin.
# The formula for this probability is P = pi / (2 * log(R)).

# In this problem, the starting point is (0,1), a neighbor of the origin.
# The target set is centered around (3600,0), so we use R = 3600.

R = 3600

# Calculate log(R)
log_R = math.log(R)

# Calculate the probability
prob = math.pi / (2 * log_R)

# We need to output the equation with the final numbers
# that leads to the result.
# The equation is: P = pi / (2 * log(3600))
# Let's print the values used.
print(f"The calculation is based on the formula: P = pi / (2 * log(R))")
print(f"Here, R = {R}")
print(f"pi is approximately {math.pi}")
print(f"log(R) = log({R}) is approximately {log_R}")
print(f"So, P = {math.pi} / (2 * {log_R})")
print(f"P = {math.pi} / {2 * log_R}")
print(f"The probability is approximately: {prob}")

# Round to two significant digits for the final answer
# This is mainly for the final boxed answer. The printout shows more precision.
prob_rounded = float(f"{prob:.2g}")
print(f"\nThe approximate answer with two significant digits is: {prob_rounded}")

# Final Answer format requires the numerical value.
# The thinking process and calculation lead to this value.
# The required format is <<<value>>>
# Let's add it at the end.
# print(f"<<<{prob_rounded}>>>")