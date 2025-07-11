# Step 1: Define the Betti numbers for the moduli space of equilateral pentagons.
# These numbers are known from established results in topology.
b0 = 1
b1 = 8
b2 = 1

print("The configuration space is a smooth surface whose topology is described by its Betti numbers.")
print(f"b_0 = {b0} (the number of connected components)")
print(f"b_1 = {b1} (related to the number of 'holes' or 'handles')")
print(f"b_2 = {b2} (the number of 'voids')")
print("-" * 30)

# Step 2: Calculate the Euler characteristic (chi) from the Betti numbers.
# The formula is chi = b0 - b1 + b2.
chi = b0 - b1 + b2

print("The Euler characteristic, chi, is calculated from the Betti numbers.")
print(f"The formula is: chi = b0 - b1 + b2")
print(f"chi = {b0} - {b1} + {b2} = {chi}")
print("-" * 30)

# Step 3: Calculate the genus (g) from the Euler characteristic.
# For a compact, orientable surface, the formula is chi = 2 - 2g.
# We can rearrange this to solve for g: g = (2 - chi) / 2.
g = (2 - chi) // 2

print("The genus, g, is found using the Euler-Poincaré formula: chi = 2 - 2g.")
print(f"Rearranging for g, we get: 2g = 2 - chi")
print(f"2g = 2 - ({chi}) = {2 - chi}")
print(f"g = {2 - chi} / 2 = {g}")
print("-" * 30)

print(f"The genus of the configuration space of the hinged pentagon is {g}.")
