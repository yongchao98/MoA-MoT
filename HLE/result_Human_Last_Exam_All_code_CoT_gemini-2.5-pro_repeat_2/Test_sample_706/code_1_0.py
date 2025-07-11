# Probabilities of edge deletion
prob_vertical_deleted = 1/2
prob_upper_horizontal_deleted = 1/3

# In the limit c -> infinity, the random walk becomes deterministic, following
# a greedy path: right > vertical > left.

# A walker on the lower level (y=0) will stay there, moving right, with speed 1.
# However, for any finite c, the walker will eventually move to the upper level.
# Once on the upper level, it is driven to the rightmost end of its connected segment.
# At this end, its fate is determined by the existence of a vertical edge to escape back
# to the lower level.

# Probability that the vertical escape edge exists
prob_escape = 1 - prob_vertical_deleted

# Probability that the walker is trapped (vertical edge is missing)
prob_trap = prob_vertical_deleted

# If the walker escapes, it returns to the lower level and its asymptotic speed is 1.
speed_if_escape = 1

# If the walker is trapped, it gets stuck in a loop and its asymptotic speed is 0.
speed_if_trap = 0

# The limiting speed v is the expected value of these two outcomes.
# The specific value of prob_upper_horizontal_deleted does not affect the final
# result, as long as it's greater than 0, which ensures the upper level
# has "ends" where the walker's fate is decided.
v_limit = speed_if_escape * prob_escape + speed_if_trap * prob_trap

print("Let p_v be the probability that a vertical edge exists.")
print(f"p_v = 1 - (probability of deletion) = 1 - {prob_vertical_deleted} = {prob_escape}")
print("\nThe asymptotic speed depends on the outcome at the end of the first upper-level segment the walker explores.")
print("With probability p_v, the walker escapes to the lower level, and the speed is 1.")
print("With probability 1-p_v, the walker is trapped, and the speed is 0.")
print("\nThe expected asymptotic speed v is therefore:")
print(f"v = ({speed_if_escape} * p_v) + ({speed_if_trap} * (1 - p_v))")
print(f"v = ({speed_if_escape} * {prob_escape}) + ({speed_if_trap} * {prob_trap}) = {v_limit}")
print("\nFinal Result:")
print(f"The limit of the asymptotic speed is {v_limit}")