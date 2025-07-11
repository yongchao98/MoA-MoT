# Plan: Determine the experimental paradigm for which each AI capability is necessary.
# The options are "both", "room", "cliff", or "neither".
# The final output will be a single string with the answers for I-V separated by hyphens.

# I. Goal-driven locomotion is needed for the "cliff" but not the "room".
capability_I = "cliff"

# II. An MST analogue for processing optic flow is the key mechanism in the "room".
capability_II = "room"

# III. Evaluating relative size is a depth cue for the "cliff".
capability_III = "cliff"

# IV. Resolving binocular disparity is a depth cue for the "cliff".
capability_IV = "cliff"

# V. Postural control is necessary to stand in the "room" and to crawl and stop in the "cliff" experiment.
capability_V = "both"

# The final answer combines each individual determination in order.
final_answer = f"{capability_I}-{capability_II}-{capability_III}-{capability_IV}-{capability_V}"

print(final_answer)