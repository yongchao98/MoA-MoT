# The task is to identify the primary parameter group that changes for four variations of playing a violin.
# The variations are: (1) sul ponticello, (2) with a bridge mute, (3) in a helium-filled room, and (4) on the E string.
# We also need to identify the direction of change for the last parameter in the group for variation (2).

# Based on the physics of the violin:
# 1. 'sul ponticello' primarily affects the string's harmonic content, which is modeled by group ii (nu).
# 2. A bridge mute primarily dampens higher-frequency body resonances, modeled by group iv (mu, a2, f2).
#    Adding mass lowers the resonance frequency f2, so the direction is 'down'.
# 3. Helium changes the speed of sound, which drastically alters the main air resonance (f1), part of group iii (a1, f1).
# 4. Changing from the A string to the E string changes the fundamental open-string frequency, F, which is group i.

# The list of groups corresponding to variations (1), (2), (3), (4)
group_identities = ["ii", "iv", "iii", "i"]

# The direction of change for the last parameter of the group for variation (2)
direction_for_2 = "down"

# The final answer is a comma-separated string of these findings.
# This code will output each part of the final answer string, joined by commas.
final_answer_parts = group_identities + [direction_for_2]

# Print the final result in the required format.
print(",".join(final_answer_parts))