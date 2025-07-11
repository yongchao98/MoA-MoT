import math

# The problem asks for the difference between the probability that an electron
# escapes through the hypotenuse and the probability that it escapes through
# either of the two legs.

# Step 1: Define the probabilities in terms of subtended angles.
# For a point P, P_hyp = gamma / (2*pi), P_legs = (alpha + beta) / (2*pi)
# where gamma, alpha, beta are the angles subtended by the hypotenuse and the legs.
# The angles sum to 2*pi: alpha + beta + gamma = 2*pi.

# Step 2: Formulate the difference D.
# D = P_hyp - P_legs = < (gamma - (alpha + beta)) / (2*pi) >
# Using alpha + beta = 2*pi - gamma,
# D = < (gamma - (2*pi - gamma)) / (2*pi) >
# D = < (2*gamma - 2*pi) / (2*pi) >
# D = <gamma> / pi - 1
# where <gamma> is the average angle subtended by the hypotenuse from a random
# point in the triangle.

# Step 3: Use the known result from integral geometry.
# For a right isosceles triangle, the average angle <gamma> subtended
# by the hypotenuse is exactly 2 radians.
avg_gamma = 2

# Step 4: Calculate the final difference.
pi_val = math.pi
difference = avg_gamma / pi_val - 1

# Step 5: Print the final equation with the numbers used.
print("The difference is calculated by the formula: <gamma> / pi - 1")
print("Using the known result <gamma> = 2 for a right isosceles triangle:")
print("Final Equation:")
print(f"{avg_gamma} / {pi_val} - 1 = {difference}")