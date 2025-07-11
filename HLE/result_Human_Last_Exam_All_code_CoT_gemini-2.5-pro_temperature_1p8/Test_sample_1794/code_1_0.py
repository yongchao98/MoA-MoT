# Part 1: How many distinct planar (non-crossing) graphs are there in the Feynman diagrams,
# excluding the diagrams with vertex corrections?

# Based on the analysis of 3-loop, 4-point, planar, primitive Feynman diagrams in phi^3 theory,
# there are two such distinct topologies. These are known as the "ladder" and "non-ladder" (or "tennis-court") graphs.
num_graphs = 2

print("--- Answers ---")
print(f"1. The number of distinct planar graphs is: {num_graphs}")
print("")

# Part 2: At d=4-2epsilon spacetime dimensions, what is the number of the power of the
# leading divergent term of the epsilon expansion of the Feynman integral of the massless
# on-shell diagram near d=4 dimension?

# For a massless on-shell L-loop diagram, the leading divergences are of infrared origin.
# Explicit calculations for these diagrams show that the leading pole in epsilon has the form 1/epsilon^(2L).
loop_order = 3

# The final equation for the power of the leading term is: Power = -2 * L
power_of_divergence = -2 * loop_order

# We print the equation and the final result as requested.
print("2. The power of the leading divergent term is calculated as:")
print(f"   Power = -2 * {loop_order} = {power_of_divergence}")