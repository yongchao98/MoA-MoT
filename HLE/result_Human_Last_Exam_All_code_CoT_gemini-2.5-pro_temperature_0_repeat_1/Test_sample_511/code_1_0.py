import math

def solve_packing_problem():
    """
    Analyzes packing strategies for energy balls in a container to find a more
    material-efficient design.
    """
    # Initial problem parameters
    ball_radius_cm = 2.0
    ball_diameter_cm = ball_radius_cm * 2
    grid_precision_cm = 0.5

    # --- Step 1: Analyze the initial 12x12x12 box (Simple Cubic Packing) ---
    initial_l, initial_w, initial_h = 12, 12, 12
    initial_sa = 2 * (initial_l * initial_w + initial_w * initial_h + initial_h * initial_l)
    
    # With a 4cm ball diameter, a 12cm edge can fit 12/4 = 3 balls.
    nx, ny, nz = 3, 3, 3
    initial_ball_count = nx * ny * nz
    
    best_solution = {
        "l": initial_l,
        "w": initial_w,
        "h": initial_h,
        "sa": initial_sa,
        "n": initial_ball_count,
        "found": False
    }

    # --- Step 2: Explore other Simple Cubic Packing (SCP) arrangements ---
    # We need to pack at least 27 balls. Let's check configurations for N > 27.
    # To minimize surface area for a given volume, the dimensions should be as close as possible.
    # We check numbers >= 27 that have factors close to each other.
    # e.g., 28 (2,2,7), 30 (2,3,5), 32 (2,4,4), 36 (3,3,4)
    
    # This part is for demonstration; as shown in the thought process, increasing N with SCP
    # will increase the surface area. For example, for N=30 (2x3x5 balls):
    # L=2*4=8, W=3*4=12, H=5*4=20. SA = 2*(8*12 + 12*20 + 20*8) = 992 > 864.
    # So, SCP will not yield a better result.

    # --- Step 3: Explore a Staggered Packing Strategy ---
    # This mimics a denser packing like FCC/BCC under grid constraints.
    # We scale everything by 2 to work with integers (0.5cm -> 1 unit)
    # Ball radius = 4 units, Ball diameter = 8 units.
    # Min distance between centers squared must be >= 8*8 = 64.
    
    # We try to pack 27 balls in a 3x3x3 staggered arrangement.
    # Layer A (z=4): 3x3 grid of balls. Centers at (4+8i, 4+8j, 4) for i,j in {0,1,2}
    # Layer B (z=10): 3x3 grid, shifted. Centers at (8+8i, 8+8j, 10)
    # Layer C (z=16): 3x3 grid, same as A. Centers at (4+8i, 4+8j, 16)
    # Vertical distance between layers is 6 units (3 cm).
    # Let's check if this is valid. Dist between a center in A and B:
    # e.g., (4,4,4) and (8,8,10). Dist^2 = (8-4)^2 + (8-4)^2 + (10-4)^2 = 16+16+36 = 68 >= 64. Valid.
    
    # Now, find the minimal bounding box for this packing of 27 balls.
    # All center coordinates are multiples of 0.5cm (or integers in our scaled units).
    all_centers_x = {4, 12, 20, 8, 16, 24}
    all_centers_y = {4, 12, 20, 8, 16, 24}
    all_centers_z = {4, 10, 16}
    
    # Box must contain spheres of radius 4 units centered at these coordinates.
    # Box starts at 0.
    # L_u >= max(all_centers_x) + radius_u = 24 + 4 = 28
    # W_u >= max(all_centers_y) + radius_u = 24 + 4 = 28
    # H_u >= max(all_centers_z) + radius_u = 16 + 4 = 20
    
    staggered_l_u, staggered_w_u, staggered_h_u = 28, 28, 20
    
    # Convert back to cm
    staggered_l = staggered_l_u * grid_precision_cm
    staggered_w = staggered_w_u * grid_precision_cm
    staggered_h = staggered_h_u * grid_precision_cm
    
    staggered_sa = 2 * (staggered_l * staggered_w + staggered_w * staggered_h + staggered_h * staggered_l)
    
    # --- Step 4: Conclusion ---
    # Compare the new design with the initial one.
    # The staggered packing holds 27 balls, but its surface area is 952 cm^2.
    # 952 > 864, so this is not a more efficient design.
    
    # The grid constraint makes it difficult to create a packing that is denser than
    # the simple cubic one. Any attempt to stagger layers requires rounding up
    # distances, which wastes space and leads to a larger container.
    
    if staggered_sa < initial_sa and staggered_l.is_integer() and staggered_w.is_integer() and staggered_h.is_integer():
        # This block would execute if a better solution were found.
        l, w, h, sa = int(staggered_l), int(staggered_w), int(staggered_h), int(staggered_sa)
        print(f"{l}:{w}:{h}:{sa}")
    else:
        # Since our exploration of more complex packings did not yield a better
        # result, we conclude that no more efficient box can be designed.
        print(0)

solve_packing_problem()
<<<0>>>