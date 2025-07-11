# Position of the torch (melt)
torch_bin = 2024
# Position of the portal (escape)
portal_bin = 2025
# Starting position
start_bin = 0

# The probability of jumping from bin n to bin n+i is (1/3)^|i|.
# The probability of eventually escaping from the starting bin `p` can be calculated as:
# p = P(jump to portal) / (P(jump to portal) + P(jump to torch))
# This is because the walk is recurrent and the starting point is far from the absorbing states.

# Probability of jumping from start_bin to portal_bin in one step
# The displacement is portal_bin - start_bin
dist_portal = portal_bin - start_bin
prob_jump_portal_num = 1
prob_jump_portal_den = 3**dist_portal

# Probability of jumping from start_bin to torch_bin in one step
# The displacement is torch_bin - start_bin
dist_torch = torch_bin - start_bin
prob_jump_torch_num = 1
prob_jump_torch_den = 3**dist_torch

# The probability of escape is prob_portal / (prob_portal + prob_torch)
# p = (1/3^2025) / (1/3^2025 + 1/3^2024)
# We can express this using a common denominator
common_den = 3**2025
num_portal = 1
num_torch = 3

# The final probability calculation
# Numerator of the final probability
final_prob_num = num_portal
# Denominator of the final probability
final_prob_den = num_portal + num_torch

# We print the calculation steps
print(f"The probability of escaping is given by the ratio of the direct jump probabilities:")
print(f"P(Escape) = P(jump to {portal_bin}) / (P(jump to {portal_bin}) + P(jump to {torch_bin}))")
print(f"P(jump to {portal_bin}) is proportional to (1/3)^{dist_portal}")
print(f"P(jump to {torch_bin}) is proportional to (1/3)^{dist_torch}")
print(f"P(Escape) = (1/3)^{dist_portal} / ((1/3)^{dist_portal} + (1/3)^{dist_torch})")
print(f"P(Escape) = 1 / (1 + (1/3)^{dist_torch} / (1/3)^{dist_portal})")
print(f"P(Escape) = 1 / (1 + 3**{dist_portal} / 3**{dist_torch})")
print(f"P(Escape) = 1 / (1 + 3**({dist_portal}-{dist_torch}))")
print(f"P(Escape) = 1 / (1 + 3**{dist_portal - dist_torch})")
print(f"P(Escape) = {final_prob_num} / {final_prob_den}")
print(f"The final probability is {final_prob_num/final_prob_den}")
