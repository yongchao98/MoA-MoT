import math

def solve_candy_packing():
    """
    This function solves the candy packing problem by first confirming the problem
    formulation and then constructing a dense packing layer by layer to find the
    maximum number of spheres.
    """
    
    # The user's problem formulation is correct. The conversion from a 12x12x11 cm box
    # and 2cm radius spheres to the given integer grid coordinates ([4, 20], [4, 20], [4, 18])
    # and the non-overlapping constraint ((xi-xj)^2 + ... >= 64) is accurate.
    # Therefore, we proceed to find the maximized value n.
    
    # We will construct a dense packing by arranging spheres in layers.

    # Box and sphere parameters in grid units (where 1 unit = 0.5 cm)
    x_range, y_range, z_range = (4, 20), (4, 20), (4, 18)
    sphere_diameter = 8  # 4 cm / 0.5 cm/unit = 8 units

    # --- Layer 1 (Pattern A) ---
    # We place the first layer at the bottom of the allowed space, z = 4.
    # The spheres are arranged in a square grid in the xy-plane.
    # The step size between centers along x and y is the sphere diameter (8 units).
    x_positions_L1 = range(x_range[0], x_range[1] + 1, sphere_diameter)
    y_positions_L1 = range(y_range[0], y_range[1] + 1, sphere_diameter)
    n_layer1 = len(x_positions_L1) * len(y_positions_L1)
    
    # --- Layer 2 (Pattern B) ---
    # This layer is staggered and placed in the hollows of Layer 1.
    # A hollow is centered horizontally between four spheres of Layer 1, so the offset is 4 units.
    # The vertical distance (dz) to the next layer is the minimum integer satisfying:
    # (offset_x)^2 + (offset_y)^2 + dz^2 >= diameter^2
    # 4^2 + 4^2 + dz^2 >= 8^2  =>  16 + 16 + dz^2 >= 64  =>  dz^2 >= 32
    # So, the minimum integer vertical separation is dz = math.ceil(math.sqrt(32)) which is 6.
    z_1 = z_range[0] # 4
    z_2 = z_1 + 6   # 10
    
    x_positions_L2 = range(x_range[0] + 4, x_range[1] + 1, sphere_diameter)
    y_positions_L2 = range(y_range[0] + 4, y_range[1] + 1, sphere_diameter)
    n_layer2 = len(x_positions_L2) * len(y_positions_L2)
    
    # --- Layer 3 (Pattern A) ---
    # This layer is placed on top of Layer 2, repeating the pattern of Layer 1.
    # The vertical distance is the same, dz = 6.
    z_3 = z_2 + 6  # 10 + 6 = 16. This is within the z_range [4, 18].
    
    # The number of spheres is the same as Layer 1.
    n_layer3 = n_layer1

    # We check if a 4th layer is possible. The next layer would be at z = 16 + 6 = 22.
    # This is outside the allowed z_range [4, 18], so only 3 layers fit.

    # The problem formulation is correct: YES.
    # Now we output the maximized value n by summing the spheres in each layer.
    
    total_n = n_layer1 + n_layer2 + n_layer3

    print("The problem formulation is correct.")
    print("To find the maximum number of eyeball candies (n), we can stack them in layers:")
    print(f"Number of candies in the first layer = {n_layer1}")
    print(f"Number of candies in the second (staggered) layer = {n_layer2}")
    print(f"Number of candies in the third layer = {n_layer3}")
    print("The final equation for the total number of candies is:")
    print(f"{n_layer1} + {n_layer2} + {n_layer3} = {total_n}")

solve_candy_packing()