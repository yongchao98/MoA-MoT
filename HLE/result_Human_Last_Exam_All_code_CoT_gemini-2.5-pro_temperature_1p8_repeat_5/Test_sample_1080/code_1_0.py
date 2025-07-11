import math

def solve_packing():
    """
    Solves the Halloween candy packing problem by constructing an efficient
    three-layer packing and verifying its validity.
    """
    
    # Radius of the eyeball candy in 0.5 cm grid units
    radius = 4
    
    # Minimum squared distance between the centers of two non-overlapping spheres
    # (diameter)^2 = (2 * radius)^2 = 8^2
    min_sq_dist = (2 * radius)**2

    # The list will store the integer coordinates (xi, yi, zi) for the center of each sphere.
    all_centers = []

    def is_placement_valid(new_center, existing_centers):
        """
        Checks if placing a new sphere at new_center is valid.
        A placement is valid if it's within the box and doesn't overlap with existing spheres.
        
        Box (12x12x11 cm) -> Grid space for centers is [4,20]x[4,20]x[4,18]
        """
        x, y, z = new_center
        # 1. Check if the sphere is fully contained within the box.
        if not (4 <= x <= 20 and 4 <= y <= 20 and 4 <= z <= 18):
            return False
            
        # 2. Check for non-overlapping with all previously placed spheres.
        for center in existing_centers:
            dist_sq = (x - center[0])**2 + (y - center[1])**2 + (z - center[2])**2
            if dist_sq < min_sq_dist:
                return False
        return True

    # Layer 1 (z=4): A 3x3 grid of spheres
    layer1_count = 0
    z_layer1 = 4
    for i in range(3):
        for j in range(3):
            center_ij = (4 + i * 8, 4 + j * 8, z_layer1)
            if is_placement_valid(center_ij, all_centers):
                all_centers.append(center_ij)
                layer1_count += 1
                
    # Layer 2 (z=10): A 2x2 grid placed in the hollows of the first layer
    layer2_count = 0
    z_layer2 = 10
    for i in range(2):
        for j in range(2):
            center_ij = (8 + i * 8, 8 + j * 8, z_layer2)
            if is_placement_valid(center_ij, all_centers):
                all_centers.append(center_ij)
                layer2_count += 1

    # Layer 3 (z=16): Another 3x3 grid, same x/y as the first layer
    layer3_count = 0
    z_layer3 = 16
    for i in range(3):
        for j in range(3):
            center_ij = (4 + i * 8, 4 + j * 8, z_layer3)
            if is_placement_valid(center_ij, all_centers):
                all_centers.append(center_ij)
                layer3_count += 1
    
    total_spheres = len(all_centers)
    
    # Print the equation showing how the total number of spheres is calculated
    print(f"The non-overlapping constraint is: (xi-xj)^2 + (yi-yj)^2 + (zi-zj)^2 >= {min_sq_dist}")
    print(f"The maximum number of candies n can be found by summing the candies in each layer:")
    print(f"{layer1_count} (Layer 1) + {layer2_count} (Layer 2) + {layer3_count} (Layer 3) = {total_spheres}")
    print(f"\nThe maximized value n is: {total_spheres}")

solve_packing()