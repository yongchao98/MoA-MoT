# The dimension of a Lie group of symmetries is the number of its generators.
# We sum the number of generators for each type of symmetry transformation
# that leaves the 3D incompressible Navier-Stokes equations invariant.

# Time translation: t' = t + s
time_translation = 1

# Space translations: x_i' = x_i + c_i for i=1,2,3
space_translations = 3

# Spatial rotations (SO(3) group): x' = Rx
rotations = 3

# Galilean boosts: x' = x + vt, u' = u + v
galilean_boosts = 3

# Scaling (dilation): t' = a^2 t, x' = a x, u' = a^-1 u
scaling_symmetry = 1

# The total dimension is the sum of these individual components.
total_dimension = (
    time_translation + 
    space_translations + 
    rotations + 
    galilean_boosts + 
    scaling_symmetry
)

print(f"The dimension is the sum of the generators for each symmetry:")
print(f"Time Translation: {time_translation}")
print(f"Space Translations: {space_translations}")
print(f"Rotations: {rotations}")
print(f"Galilean Boosts: {galilean_boosts}")
print(f"Scaling Symmetry: {scaling_symmetry}")
print("-" * 30)
print(f"Total Dimension = {time_translation} + {space_translations} + {rotations} + {galilean_boosts} + {scaling_symmetry} = {total_dimension}")