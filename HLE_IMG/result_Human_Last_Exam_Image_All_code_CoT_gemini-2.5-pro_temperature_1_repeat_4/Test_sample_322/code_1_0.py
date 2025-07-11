# This script prints the final sequence based on the visual analysis of vortex trajectories.
# The analysis identifies the unique vortex in each of the 16 plots.
# An uppercase letter (R, G, B) indicates the vortex is twice as strong as the other two.
# This is identified by its smaller, more central trajectory.
# A lowercase letter (r, g, b) indicates the vortex is half as strong as the other two.
# This is identified by its larger, more extensive, or chaotic trajectory.

# The sequence is determined by analyzing each plot from 1 to 16.
# 1: G, 2: B, 3: r, 4: r, 5: B, 6: G, 7: g, 8: g
# 9: G, 10: R, 11: b, 12: b, 13: B, 14: R, 15: r, 16: g
final_sequence = "GBrrBGggGRbBbrg"

print(final_sequence)