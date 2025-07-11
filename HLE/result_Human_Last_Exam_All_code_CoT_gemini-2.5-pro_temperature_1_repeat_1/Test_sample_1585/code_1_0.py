# Step 1: Define the values for n_2 and n_3 based on the problem analysis.
# The graph satisfying the conditions is the Line Graph of the Petersen Graph.
# The number of vertices in this graph is 15.
# n_g = 15

# Step 2: Determine n_2.
# A graph is 2-planar if it has a thickness of at most 2.
# The thickness of the Line Graph of the Petersen Graph is 2.
# Assuming this is the smallest such graph, n_2 is 15.
n_2 = 15

# Step 3: Determine n_3.
# A graph is 3-planar if it has a thickness of at most 3.
# A graph with thickness 2 also has thickness 3.
# The smallest n for such a graph is therefore also 15.
n_3 = 15

# Step 4: Calculate the final result using the given formula.
result = (n_2 + n_3) * n_2

# Step 5: Print the final equation with the numbers substituted.
print(f"({n_2} + {n_3}) * {n_2} = {result}")

# Step 6: Output the final answer in the specified format.
print(f"<<<{result}>>>")