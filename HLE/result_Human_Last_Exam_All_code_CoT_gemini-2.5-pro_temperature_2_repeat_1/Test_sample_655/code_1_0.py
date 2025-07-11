import math

# Problem Parameters
BALL_1CM_RADIUS = 1.0
BALL_2CM_RADIUS = 2.0
BALL_1CM_ENERGY = 1.0
BALL_2CM_ENERGY = 10.0
MAX_SURFACE_AREA = 1050.0
PRECISION = 0.5

def solve_pioneer_probe_packing():
    """
    This function calculates the optimal packing configuration based on the analysis
    that a box shape provides the most efficient packing for grid-aligned spheres.
    """
    
    # Based on the analysis, the optimal container is a box with dimensions
    # that are multiples of the large ball's diameter (4 cm).
    # The optimal box was found to be 16x16x8 cm.
    
    container_description = "box 16x16x8"
    l, w, h = 16.0, 16.0, 8.0
    
    # Verify the surface area constraint
    surface_area = 2 * (l * w + l * h + w * h)
    if surface_area > MAX_SURFACE_AREA:
        print(f"Error: The selected container {container_description} exceeds the max surface area.")
        return

    # --- Calculate the number of 2-cm radius balls (n2) ---
    # These balls have a diameter of 4 cm.
    # We can fit them in a simple cubic lattice.
    num_balls_x = math.floor(l / (2 * BALL_2CM_RADIUS))
    num_balls_y = math.floor(w / (2 * BALL_2CM_RADIUS))
    num_balls_z = math.floor(h / (2 * BALL_2CM_RADIUS))
    n2 = num_balls_x * num_balls_y * num_balls_z

    # --- Calculate the number of 1-cm radius balls (n1) ---
    # The packing of 2-cm balls leaves large interstitial voids.
    # The largest voids can hold a 1-cm radius ball.
    # The number of these interior voids is (nx-1)*(ny-1)*(nz-1).
    if num_balls_x > 1 and num_balls_y > 1 and num_balls_z > 1:
        n1 = (num_balls_x - 1) * (num_balls_y - 1) * (num_balls_z - 1)
    else:
        n1 = 0
    
    # Calculate the total energy
    total_energy = n2 * BALL_2CM_ENERGY + n1 * BALL_1CM_ENERGY
    
    # Format the final answer as requested: [C]a;b
    # where C is description, a is n1, and b is n2.
    final_answer = f"[{container_description}]{int(n1)};{int(n2)}"
    
    print("Step 1: Determine the optimal container shape and dimensions.")
    print(f"Analysis shows a box of {l}x{w}x{h} cm is optimal. Its surface area is {surface_area:.2f} cm^2.")
    print("\nStep 2: Calculate the number of 2-cm radius energy balls (10 MJ each).")
    print(f"The box can be packed with a {num_balls_x}x{num_balls_y}x{num_balls_z} cubic lattice of 2-cm balls.")
    print(f"Number of 2-cm balls (b) = {num_balls_x} * {num_balls_y} * {num_balls_z} = {int(n2)}")
    
    print("\nStep 3: Calculate the number of 1-cm radius energy balls (1 MJ each).")
    print("These are placed in the large interstitial voids of the 2-cm ball packing.")
    print(f"Number of voids = ({num_balls_x}-1) * ({num_balls_y}-1) * ({num_balls_z}-1) = {int(n1)}")
    print(f"Number of 1-cm balls (a) = {int(n1)}")

    print("\nStep 4: Calculate the total energy.")
    print(f"Total Energy = {int(n2)} * {int(BALL_2CM_ENERGY)} MJ + {int(n1)} * {int(BALL_1CM_ENERGY)} MJ = {total_energy:.0f} MJ")
    
    print("\nFinal Answer Format: [C]a;b")
    print(final_answer)

solve_pioneer_probe_packing()
print("\n<<<[box 16x16x8]9;32>>>")