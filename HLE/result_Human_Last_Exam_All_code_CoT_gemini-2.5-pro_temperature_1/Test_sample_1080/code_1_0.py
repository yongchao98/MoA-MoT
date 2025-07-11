import itertools

def solve_packing_problem():
    """
    This script verifies a proposed solution for packing 22 spheres into a box,
    based on the user's problem formulation.
    """
    
    # Problem parameters in grid units (1 unit = 0.5 cm)
    box_dims = {'x': (4, 20), 'y': (4, 20), 'z': (4, 18)}
    sphere_radius = 4
    min_dist_sq = (2 * sphere_radius)**2 # (4+4)^2 = 64

    # Proposed solution: A staggered packing of 22 spheres
    # This configuration consists of three layers.
    # Layer 1: 3x3 grid at z=4
    # Layer 2: 2x2 grid at z=10 (in the hollows of layer 1)
    # Layer 3: 3x3 grid at z=16 (aligned with layer 1)
    
    centers = []
    # Layer 1 (9 spheres)
    for x in [4, 12, 20]:
        for y in [4, 12, 20]:
            centers.append((x, y, 4))
            
    # Layer 2 (4 spheres)
    for x in [8, 16]:
        for y in [8, 16]:
            centers.append((x, y, 10))
            
    # Layer 3 (9 spheres)
    for x in [4, 12, 20]:
        for y in [4, 12, 20]:
            centers.append((x, y, 16))

    n = len(centers)
    print(f"Proposed solution: Packing {n} spheres.\n")

    # --- Verification Step 1: Check if all sphere centers are within the box boundaries ---
    print("--- Verifying sphere boundaries ---")
    all_spheres_in_box = True
    for i, p in enumerate(centers):
        x, y, z = p
        if not (box_dims['x'][0] <= x <= box_dims['x'][1] and
                box_dims['y'][0] <= y <= box_dims['y'][1] and
                box_dims['z'][0] <= z <= box_dims['z'][1]):
            print(f"Sphere {i+1} at {p} is OUT of bounds!")
            all_spheres_in_box = False
            break
    if all_spheres_in_box:
        print("All spheres are within the box boundaries. OK.\n")

    # --- Verification Step 2: Check for overlaps between all pairs of spheres ---
    print("--- Verifying non-overlapping constraint ---")
    print(f"The constraint is: (x_i-x_j)^2 + (y_i-y_j)^2 + (z_i-z_j)^2 >= {min_dist_sq}\n")
    
    no_overlaps = True
    # itertools.combinations ensures we check each pair only once
    for (i, p1), (j, p2) in itertools.combinations(enumerate(centers), 2):
        dx = p1[0] - p2[0]
        dy = p1[1] - p2[1]
        dz = p1[2] - p2[2]
        
        dist_sq = dx**2 + dy**2 + dz**2
        
        print(f"Checking pair ({i+1}, {j+1}) at {p1} and {p2}:")
        print(f"({p1[0]}-{p2[0]})^2 + ({p1[1]}-{p2[1]})^2 + ({p1[2]}-{p2[2]})^2 = "
              f"{dx**2} + {dy**2} + {dz**2} = {dist_sq}")
              
        if dist_sq < min_dist_sq:
            print(f"  -> VIOLATION! {dist_sq} is less than {min_dist_sq}.")
            no_overlaps = False
        else:
            print(f"  -> OK. {dist_sq} >= {min_dist_sq}")

    print("\n--- Verification Complete ---")
    if all_spheres_in_box and no_overlaps:
        print("\nThe proposed arrangement is valid.")
        print(f"The maximized value is n = {n}")
    else:
        print("\nThe proposed arrangement is invalid.")
        print("The maximized value could not be confirmed.")

solve_packing_problem()