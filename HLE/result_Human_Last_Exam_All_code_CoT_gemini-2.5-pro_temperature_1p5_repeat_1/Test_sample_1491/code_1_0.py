# The analysis for each capability is as follows:
# I. Goal-driven locomotion triggered by visual stimuli -> both
# II. A functional analogue of cortical area MST in humans -> room
# III. A perceptual mechanism for evaluating the relative size of elements -> cliff
# IV. A mechanism for resolving binocular disparity -> cliff
# V. A mechanism for controlling the effectors that maintain posture -> room

# This script concatenates the answers into the required hyphen-separated format.
answer_i = "both"
answer_ii = "room"
answer_iii = "cliff"
answer_iv = "cliff"
answer_v = "room"

final_answer = "-".join([answer_i, answer_ii, answer_iii, answer_iv, answer_v])

print(final_answer)