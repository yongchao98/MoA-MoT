import math

def solve_eyeball_packing():
    """
    This function solves the sphere packing problem by constructing a highly efficient
    layered packing arrangement and then summing the number of spheres.
    """
    
    # Box and candy properties in grid units (1 unit = 0.5 cm)
    # Box: 24x24x22 units. Radius: 4 units. Diameter: 8 units.
    
    # Layer 1: A 3x3 grid at the bottom of the box.
    # The lowest possible z-coordinate for a center is 4 (real coordinate 2 cm).
    layer1_coords = []
    z1 = 4
    # The centers are placed 8 units (4 cm) apart.
    # We choose coordinates {4, 12, 20} which correspond to {2, 6, 10} cm.
    # These spheres will fit within the 12x12 cm (24x24 units) dimensions.
    for x in [4, 12, 20]:
        for y in [4, 12, 20]:
            layer1_coords.append((x, y, z1))
    n1 = len(layer1_coords)

    # Layer 2: A 2x2 grid placed in the hollows of Layer 1.
    # The non-overlapping constraint is (dx^2 + dy^2 + dz^2 >= 8^2).
    # For a sphere in the hollow, min dx=4, min dy=4. So 4^2 + 4^2 + dz^2 >= 64.
    # This means 32 + dz^2 >= 64 -> dz^2 >= 32 -> dz >= sqrt(32) approx 5.66.
    # Since z must be an integer, the minimum difference in z-coordinates is 6.
    layer2_coords = []
    z2 = z1 + 6  # 4 + 6 = 10
    # The centers for the hollows are at {8, 16}.
    for x in [8, 16]:
        for y in [8, 16]:
            layer2_coords.append((x, y, z2))
    n2 = len(layer2_coords)
    
    # Layer 3: Another 3x3 grid, stacked on top of Layer 2.
    # The z-separation logic is the same as between layers 1 and 2.
    layer3_coords = []
    z3 = z2 + 6 # 10 + 6 = 16
    # This layer has the same xy-coordinates as the first layer.
    # We check if z3=16 is valid. Max z is 18, so it is.
    for x in [4, 12, 20]:
        for y in [4, 12, 20]:
            layer3_coords.append((x, y, z3))
    n3 = len(layer3_coords)
    
    total_spheres = n1 + n2 + n3

    print(f"A dense packing can be achieved by stacking three layers:")
    print(f"Layer 1 is a 3x3 grid and contains {n1} spheres.")
    print(f"Layer 2 is a 2x2 grid and contains {n2} spheres.")
    print(f"Layer 3 is a 3x3 grid and contains {n3} spheres.")
    print("\nThe final equation for the total number of spheres is:")
    print(f"{n1} + {n2} + {n3} = {total_spheres}")

# Run the solver
solve_eyeball_packing()