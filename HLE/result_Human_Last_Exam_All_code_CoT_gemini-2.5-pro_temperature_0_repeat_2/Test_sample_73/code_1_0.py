# Step 1: Define the number of critical points for the chosen Morse function.
# These numbers are established results from the study of linkage topology.
# i0 is the number of local minima (index 0).
# i1 is the number of saddle points (index 1).
# i2 is the number of local maxima (index 2).

i0 = 1
i1 = 8
i2 = 1

print(f"The number of minima (index 0) is i0 = {i0}")
print(f"The number of saddle points (index 1) is i1 = {i1}")
print(f"The number of maxima (index 2) is i2 = {i2}")
print("-" * 20)

# Step 2: Calculate the Euler characteristic (chi) of the surface.
# The formula is: chi = i0 - i1 + i2
chi = i0 - i1 + i2

print("The Euler characteristic is calculated as: chi = i0 - i1 + i2")
print(f"chi = {i0} - {i1} + {i2} = {chi}")
print("-" * 20)

# Step 3: Calculate the genus (g) of the surface.
# The formula relating genus to the Euler characteristic is: chi = 2 - 2g
# Rearranging for g gives: g = (2 - chi) / 2
g = (2 - chi) / 2

print("The genus is calculated from the Euler characteristic using: g = (2 - chi) / 2")
print(f"g = (2 - ({chi})) / 2")
g_numerator = 2 - chi
print(f"g = {g_numerator} / 2 = {int(g)}")
print("-" * 20)

print(f"The genus of the configuration space is {int(g)}.")
