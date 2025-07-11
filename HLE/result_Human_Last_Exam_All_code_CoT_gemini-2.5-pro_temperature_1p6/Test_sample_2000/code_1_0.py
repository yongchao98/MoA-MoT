import math

# The number of hyperedges in the hypergraph.
m = 3

# The formula for the maximum generalized hypertreewidth (ghw) of a hypergraph
# with m hyperedges is floor(m/2) + 1.
floor_val = math.floor(m / 2)
result = floor_val + 1

# Print the final equation with all its components.
print(f"The number of hyperedges (m) is: {m}")
print(f"The formula for the maximum ghw is: floor(m / 2) + 1")
print(f"Step 1: Calculate m / 2 = {m/2}")
print(f"Step 2: Calculate floor(m / 2) = floor({m/2}) = {floor_val}")
print(f"Step 3: Calculate floor(m / 2) + 1 = {floor_val} + 1 = {result}")
print(f"The maximum generalised hypertreewidth is: {result}")