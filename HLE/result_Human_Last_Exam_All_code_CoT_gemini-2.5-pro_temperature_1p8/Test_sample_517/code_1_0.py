# The problem asks for the limit of the probability p_n as n approaches infinity.
# This can be solved using symmetry arguments about the asymptotic direction of the walk.

# 1. Probability of the walk escaping into the right half-plane {x > 0}.
# Due to the symmetry of the walk's rules with respect to the y-axis,
# there is an equal chance of escaping to the left or right.
prob_half_plane = 1 / 2

# 2. Number of quadrants in the right half-plane.
# The probability of escaping into the right half-plane is split equally
# between the upper-right quadrant {x > 0, y > 0} and the
# lower-right quadrant {x > 0, y < 0}.
num_quadrants_in_half_plane = 2

# 3. The limit p_n is the probability of the walk escaping along the positive x-axis,
# which corresponds to the share of one of these two quadrants.
final_prob = prob_half_plane / num_quadrants_in_half_plane

# Output the components of the final calculation.
print("The calculation is based on the following reasoning:")
print(f"The probability that the walk eventually enters the right half-plane is {prob_half_plane}.")
print(f"This probability is distributed over {num_quadrants_in_half_plane} quadrants within that half-plane.")
print(f"The final limiting probability is the ratio: {prob_half_plane} / {num_quadrants_in_half_plane} = {final_prob}")