import math

# --- Problem Parameters ---
sphere_radius = 2.0
cylinder_radius = 10.0
cylinder_height = 11.0

# --- Heuristic Calculation ---

# 1. Calculate spheres per layer (based on a hexagonal packing structure)
# This is a known result for this ratio of radii (R_eff/r = 4)
spheres_per_layer = 19 

# 2. Calculate the number of layers
sphere_diameter = 2 * sphere_radius
# Vertical distance between centers in a dense packing
layer_separation = sphere_diameter * math.sqrt(2.0/3.0)
# The available height for sphere centers
available_height = cylinder_height - 2 * sphere_radius
# Number of layers = floor(available_height / separation) + 1
# which gives 3 layers as determined in the explanation
num_layers = 3

# 3. Calculate total number of spheres
total_spheres = spheres_per_layer * num_layers

print("This calculation is based on a geometric heuristic assuming a layered packing.")
print(f"Number of layers that can be stacked: {num_layers}")
print(f"Number of spheres per layer: {spheres_per_layer}")
print("\nFinal calculation:")
print(f"{num_layers} layers * {spheres_per_layer} spheres/layer = {total_spheres} spheres")
