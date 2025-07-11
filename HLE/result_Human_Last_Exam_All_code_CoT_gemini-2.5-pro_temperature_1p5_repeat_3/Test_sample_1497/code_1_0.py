import sys

# Based on the analysis of the violin waveform equation and the physical
# effects of each variation, we determine the corresponding parameter group.

# (1) 'sul ponticello' (bowing near the bridge) emphasizes high harmonics,
#     which is controlled by the bowing parameter nu. This maps to group ii.
group_1 = "ii"

# (2) A bridge mute adds mass to the bridge, which primarily affects the
#     body's higher-frequency structural resonances. This maps to group iv.
group_2 = "iv"

# (3) Helium in the violin body increases the speed of sound, which raises
#     the frequency of the main air resonance (f_1). This maps to group iii.
group_3 = "iii"

# (4) Changing from the A string to the E string changes the fundamental
#     frequency of the open string (F). This maps to group i.
group_4 = "i"

# For variation (2) (bridge mute), adding mass to the bridge system lowers
# the frequency of its resonances. The last parameter in group iv is f_2,
# which represents a resonance frequency. Therefore, its direction of change is down.
direction_for_2 = "down"

# The final answer is the comma-separated sequence of the groups in order,
# followed by the direction for variation (2).
final_answer = ",".join([group_1, group_2, group_3, group_4, direction_for_2])

print(final_answer)
