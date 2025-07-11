# Step 1: Define the number of critical points based on established mathematical results.
# These correspond to specific configurations of the hinged pentagon.

# Number of maxima of the area function (regular convex pentagons)
num_maxima = 2

# Number of minima of the area function (regular star-shaped pentagons)
num_minima = 2

# Number of saddle points of the area function (trapezoidal configurations)
num_saddles = 8

# Step 2: Calculate the Euler characteristic (chi) of the surface.
# The formula is chi = (#maxima + #minima) - #saddles.
chi = num_maxima + num_minima - num_saddles

print(f"The number of maxima is: {num_maxima}")
print(f"The number of minima is: {num_minima}")
print(f"The number of saddles is: {num_saddles}")
print(f"The Euler characteristic chi is calculated as {num_maxima} + {num_minima} - {num_saddles} = {chi}")

# Step 3: Calculate the genus (g) of the surface.
# The formula relating genus to the Euler characteristic for a closed orientable surface is chi = 2 - 2g.
# We can rearrange this to solve for g: g = (2 - chi) / 2.
genus = (2 - chi) / 2

print(f"The genus g is calculated from chi using the formula g = (2 - chi) / 2.")
print(f"So, g = (2 - ({chi})) / 2 = {int(genus)}")

print(f"\nThe genus of the configuration space is {int(genus)}.")
