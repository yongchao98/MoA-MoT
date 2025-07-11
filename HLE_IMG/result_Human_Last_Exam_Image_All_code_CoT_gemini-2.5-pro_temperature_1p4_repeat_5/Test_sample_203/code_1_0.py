import math

# Based on the analysis, the polygon cannot be surrounded by even a single layer
# of copies of itself because no two copies can be placed edge-to-edge without gaps.
# The image itself demonstrates this with the gaps between the vertically arranged polygons.
# When the first layer cannot be formed, the Heesch number is 0.
# The three polygons in the figure are identical.

heesch_number_1 = 0
heesch_number_2 = 0
heesch_number_3 = 0

# The problem asks for the answer as a comma-separated list for the three polygons.
print(f"{heesch_number_1}, {heesch_number_2}, {heesch_number_3}")