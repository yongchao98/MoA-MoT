# Number of points of each color
n_R = 4  # Red points
n_G = 2  # Green points
n_Y = 2  # Yellow points

# Calculate the total number of points n
n = n_R + n_G + n_Y

# This configuration is valid:
# 1. With 4 red points arranged as a convex quadrilateral, 2 green points
#    are sufficient to be inside all 4 red triangles.
# 2. With only 2 green points, no green triangles can be formed,
#    so the second condition is vacuously true.
# 3. With only 2 yellow points, no yellow triangles can be formed,
#    so the third condition is vacuously true.

# We will print the equation for the total number of points.
print(f"{n_R} + {n_G} + {n_Y} = {n}")