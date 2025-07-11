# This script prints the final 16-letter sequence identifying the unique vortex in each plot.
# The logic for identifying each vortex is as follows:
# - Uppercase (R, G, B): The vortex is twice as strong. Its trajectory is typically simpler, more central, or acts as a pivot.
# - Lowercase (r, g, b): The vortex is half as strong. Its trajectory is typically more complex, chaotic, or appears "dragged along" by the other two.

# Sequence construction based on visual analysis of the 16 plots:
# Plot 1: G (Green is central pivot)
# Plot 2: B (Blue is central pivot)
# Plot 3: R (Red is the chaotic anchor, more constrained than others)
# Plot 4: R (Red's path is simplest)
# Plot 5: B (Blue is central pivot)
# Plot 6: R (Red is central pivot)
# Plot 7: r (Red's path is most chaotic)
# Plot 8: b (Blue's path is complex, not a pivot)
# Plot 9: R (Red is central pivot)
# Plot 10: G (Green is central pivot)
# Plot 11: b (Blue's path is most chaotic/space-filling)
# Plot 12: g (Green's path is complex, not a pivot)
# Plot 13: b (Blue is entrained by the stronger R-G pair)
# Plot 14: b (Blue is entrained by the stronger R-G pair)
# Plot 15: g (Green's path is most chaotic)
# Plot 16: r (Red's path is complex, not a pivot)

final_sequence = "GBRRBRrbRGbgbbgr"
print(final_sequence)