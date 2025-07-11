import sys
import io

# Buffer to capture the output
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

# This script is designed to output the solution to the vortex trajectory analysis problem.
# The analysis is performed visually based on the principles of vortex dynamics.

# Principles of analysis:
# 1. Identify the vortex with the "odd one out" trajectory.
# 2. If the unique vortex's path is dominant, smooth, and acts as a center of rotation
#    for the other two, it is considered to have twice the strength (Uppercase: R, G, B).
# 3. If the unique vortex's path is submissive, confined, or chaotically tossed around by
#    the other two (who often form a symmetric pair), it is considered to have half the
#    strength (Lowercase: r, g, b).

# Analysis results for each plot from 1 to 16:
# 1: Green is confined by the Red/Blue pair -> weak (g)
# 2: Red has a large, simple path that the others orbit -> dominant (R)
# 3: Green/Blue form a chaotic pair, tossing the Red one around -> weak (r)
# 4: Blue has the largest, most stable orbit, enclosing the others -> dominant (B)
# 5: Red/Green form a pair that contains the Blue vortex -> weak (b)
# 6: Green is confined in a central spiral by the Red/Blue pair -> weak (g)
# 7: Red/Blue form a chaotic pair that tosses Green around -> weak (g)
# 8: Red's path is dominant, with Green/Blue orbiting it -> dominant (R)
# 9: Green has a small, weak path inside the Red/Blue pair -> weak (g)
# 10: Blue is a dominant central vortex that Red/Green clearly orbit -> dominant (B)
# 11: Red/Blue form a chaotic pair tossing Green between them -> weak (g)
# 12: Red's path is dominant, with the G/B pair orbiting inside -> dominant (R)
# 13: Blue has a very small, weak trajectory inside the Red/Green pair -> weak (b)
# 14: Similar to 13, Blue is weak and confined by the Red/Green pair -> weak (b)
# 15: Green/Blue form a chaotic pair tossing Red between them -> weak (r)
# 16: Similar to 10, Blue is a dominant central vortex orbited by Red/Green -> dominant (B)

# The final sequence is constructed by concatenating the individual results.
final_sequence = "gRrBbggRgBgRbbrB"

print(final_sequence)

# Restore stdout and get the captured output
sys.stdout = old_stdout
output = captured_output.getvalue()

# Final answer format
print(output, end="")
print(f"<<<{output.strip()}>>>")