import math

def solve():
    """
    Solves the energy ball packing problem by finding a more efficient container.
    """
    # Step 1: Analyze the initial container
    initial_side = 12.0
    initial_surface_area = 6 * initial_side**2
    
    # In a simple cubic packing, balls of diameter 4cm are placed on a grid.
    balls_per_side = math.floor(initial_side / 4.0)
    initial_ball_capacity = int(balls_per_side**3)

    print("--- Initial Container Analysis ---")
    print(f"The initial container is a {initial_side}x{initial_side}x{initial_side} cm cube.")
    print(f"The surface area is 6 * {initial_side}^2 = {initial_surface_area} cm^2.")
    print(f"With a simple grid packing, it can hold {initial_ball_capacity} balls (d=4 cm).")
    print("-" * 34)

    # Step 2 & 3: Design a new packing and container
    # We can pack 27 balls more densely using a shifted 3-layer structure (A-B-A).
    # Each layer is a 3x3 grid of 9 balls.
    # The vertical distance between layer centers (dz) and horizontal shifts (s, t)
    # must satisfy: s^2 + t^2 + dz^2 >= 4^2.
    # To minimize the overall container height, we make the vertical distance between
    # the first and third layer exactly 4cm. This implies the center layer is 2cm
    # above the first, so dz=2cm.
    # This leads to the condition for the horizontal shifts: s^2 + t^2 >= 12.
    # We choose s=3.5 and t=0.5 (multiples of 0.5) to minimize the final surface area.
    
    ball_radius = 2.0
    s_shift = 3.5
    t_shift = 0.5
    
    # Layer centers definition
    grid_coords = [2.0, 6.0, 10.0]
    
    # Layer 1 at z=2.0
    layer1_centers = [(x, y, 2.0) for x in grid_coords for y in grid_coords]
    
    # Layer 2 at z=4.0, shifted by s and t
    layer2_centers = [(x + s_shift, y + t_shift, 4.0) for x in grid_coords for y in grid_coords]
    
    # Layer 3 at z=6.0, same horizontal grid as Layer 1
    layer3_centers = [(x, y, 6.0) for x in grid_coords for y in grid_coords]

    all_centers = layer1_centers + layer2_centers + layer3_centers
    new_ball_capacity = len(all_centers)

    # Step 4: Calculate the bounding box for the new packing
    min_x = min(c[0] for c in all_centers) - ball_radius
    max_x = max(c[0] for c in all_centers) + ball_radius
    min_y = min(c[1] for c in all_centers) - ball_radius
    max_y = max(c[1] for c in all_centers) + ball_radius
    min_z = min(c[2] for c in all_centers) - ball_radius
    max_z = max(c[2] for c in all_centers) + ball_radius

    new_L = max_x - min_x
    new_W = max_y - min_y
    new_H = max_z - min_z

    # Step 5: Calculate the new surface area
    new_surface_area = 2 * (new_L * new_W + new_L * new_H + new_W * new_H)

    print("\n--- Proposed New Container ---")
    print(f"A new packing configuration for {new_ball_capacity} balls has been designed.")
    print(f"The bounding box required for this packing is {new_L}x{new_W}x{new_H} cm.")
    
    print("\n--- Final Answer Calculation ---")
    print(f"New Surface Area = 2 * (L*W + L*H + W*H)")
    print(f"New Surface Area = 2 * ({new_L}*{new_W} + {new_L}*{new_H} + {new_W}*{new_H})")
    print(f"New Surface Area = 2 * ({new_L*new_W} + {new_L*new_H} + {new_W*new_H})")
    print(f"New Surface Area = 2 * ({new_L*new_W + new_L*new_H + new_W*new_H})")
    print(f"New Surface Area = {new_surface_area}")

    if new_surface_area < initial_surface_area:
        print(f"\nSince {new_surface_area} < {initial_surface_area}, this is a more efficient design.")
        answer_d = new_surface_area
        answer_X = f"box {new_L}x{new_W}x{new_H}"
        final_answer = f"{answer_d}[{answer_X}]"
    else:
        final_answer = "0"
        
    print(f"\n<<<{final_answer}>>>")

solve()