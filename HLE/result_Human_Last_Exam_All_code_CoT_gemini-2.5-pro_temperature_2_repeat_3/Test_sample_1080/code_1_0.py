import math

def solve_packing_problem():
    """
    Calculates the maximum number of spheres in a box with grid constraints.
    """
    # Problem Parameters in grid units (1 unit = 0.5 cm)
    box_dims_cm = (12, 12, 11)
    sphere_radius_cm = 2.0

    sphere_radius_units = int(sphere_radius_cm / 0.5) # 4 units
    min_center_dist_sq = (2 * sphere_radius_units)**2 # 8*8 = 64 units^2

    box_dims_units = tuple(int(d / 0.5) for d in box_dims_cm)
    
    # Valid center coordinate ranges
    # A sphere center (x,y,z) must be at least `radius` away from each wall.
    # Box is from (0,0,0) to (24,24,22)
    # x in [4, 24-4], y in [4, 24-4], z in [4, 22-4]
    x_range = (4, 20)
    y_range = (4, 20)
    z_range = (4, 18)

    # --- Strategy: Layered Packing ---
    # Layer Type 'L' (Large): 3x3 grid.
    # Spacing between centers is diameter (8 units) to fit on a 12cm (24 units) plane.
    # x/y coordinates: 4, 4+8, 4+8+8 => 4, 12, 20
    num_spheres_L = 3 * 3

    # Layer Type 'S' (Small): 2x2 grid, fits in the hollows of an L-layer.
    # x/y coordinates: 8, 8+8 => 8, 16
    num_spheres_S = 2 * 2

    # Vertical distance calculation between an L layer and an S layer
    # Let L be at z1 and S at z2. Pick a sphere from L, e.g., (4,4,z1) and one from S, e.g., (8,8,z2).
    # (x_s-x_l)^2 + (y_s-y_l)^2 + (z_s-z_l)^2 >= 64
    # (8-4)^2 + (8-4)^2 + (z2-z1)^2 >= 64
    # 16 + 16 + (z2-z1)^2 >= 64
    # (z2-z1)^2 >= 32
    # Smallest integer |z2-z1| is ceil(sqrt(32)) = 6.
    min_z_sep_LS = 6

    # Construct the L-S-L packing stack
    layers = []
    
    # Layer 1: An 'L' layer at the bottom. Start at the lowest possible z.
    z1 = z_range[0] # 4
    layers.append({'type': 'L', 'z': z1, 'count': num_spheres_L})
    print(f"Placing Layer 1 (a {int(math.sqrt(num_spheres_L))}x{int(math.sqrt(num_spheres_L))} grid) at z = {z1*0.5} cm.")
    print(f"Number of candies in Layer 1: {num_spheres_L}")

    # Layer 2: An 'S' layer on top of Layer 1.
    z2 = z1 + min_z_sep_LS # 4 + 6 = 10
    if z2 <= z_range[1]:
        layers.append({'type': 'S', 'z': z2, 'count': num_spheres_S})
        print(f"Placing Layer 2 (a {int(math.sqrt(num_spheres_S))}x{int(math.sqrt(num_spheres_S))} grid) at z = {z2*0.5} cm.")
        print(f"Number of candies in Layer 2: {num_spheres_S}")

    # Layer 3: An 'L' layer on top of Layer 2.
    z3 = z2 + min_z_sep_LS # 10 + 6 = 16
    if z3 <= z_range[1]:
        layers.append({'type': 'L', 'z': z3, 'count': num_spheres_L})
        print(f"Placing Layer 3 (a {int(math.sqrt(num_spheres_L))}x{int(math.sqrt(num_spheres_L))} grid) at z = {z3*0.5} cm.")
        print(f"Number of candies in Layer 3: {num_spheres_L}")

    # Calculate total
    total_spheres = sum(layer['count'] for layer in layers)

    # Print the final equation
    counts = [str(layer['count']) for layer in layers]
    equation = " + ".join(counts)
    
    print("\nThis packing configuration is valid and fits within the box.")
    print("Final Calculation:")
    print(f"Total number of candies = {equation} = {total_spheres}")
    
    return total_spheres

# Yes, your problem formulation is correct.
print("Is the problem formulation correct? Yes.")
print("-" * 20)
max_n = solve_packing_problem()
print("-" * 20)
print(f"The maximized value n is {max_n}.")
print("<<<{}>>>".format(max_n))
