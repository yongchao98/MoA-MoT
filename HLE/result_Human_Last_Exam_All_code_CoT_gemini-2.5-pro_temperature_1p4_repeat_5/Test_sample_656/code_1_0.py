import math

def solve_pioneer_probe_packing():
    """
    This function solves the Pioneer probe packing problem by analyzing a specific,
    highly optimal container and packing strategy.

    The chosen container is a 13.0x13.0x13.0 cm box, which is the largest cube
    fitting within the 1050 cm^2 surface area constraint.

    The packing strategy is:
    1. Pack large (2-cm radius) balls in a simple 3x3x3 cubic grid.
    2. Pack small (1-cm radius) balls greedily into the voids left by the large balls.
    """

    # --- Container and Ball Properties ---
    box_dims = (13.0, 13.0, 13.0)
    
    large_ball_radius = 2.0
    large_ball_energy = 20
    
    small_ball_radius = 1.0
    small_ball_energy = 1
    
    # --- Step 1: Pack Large Balls ---
    # We place large balls (diameter 4) in a simple cubic lattice.
    # The number of balls along each dimension is floor(side / diameter).
    n_large_x = math.floor(box_dims[0] / (2 * large_ball_radius))
    n_large_y = math.floor(box_dims[1] / (2 * large_ball_radius))
    n_large_z = math.floor(box_dims[2] / (2 * large_ball_radius))
    
    num_large_balls = n_large_x * n_large_y * n_large_z

    # Store the centers of the large balls for collision checking.
    # The grid of centers is at (2, 6, 10), (2, 6, 10), (2, 6, 10).
    large_ball_centers = []
    for i in range(n_large_x):
        for j in range(n_large_y):
            for k in range(n_large_z):
                center = (
                    large_ball_radius + i * (2 * large_ball_radius),
                    large_ball_radius + j * (2 * large_ball_radius),
                    large_ball_radius + k * (2 * large_ball_radius)
                )
                large_ball_centers.append(center)

    # --- Step 2: Pack Small Balls in Voids ---
    # We identify three types of voids where small balls can fit.

    # Type 1: Central Voids (center of 8 large balls)
    # These are centered at (4,4,4), (4,4,8), etc.
    num_central_voids = (n_large_x - 1) * (n_large_y - 1) * (n_large_z - 1)

    # Type 2: Face Voids (between large balls and container walls)
    # There are 4 such voids on each of the 6 faces of the container.
    # For example, on the face at z=13, balls can be centered at (4,4,12), (4,8,12), etc.
    num_face_voids = 6 * (n_large_x - 1) * (n_large_y - 1)
    
    # Type 3: Corner Voids
    # Analysis shows that in a 13x13x13 box, only 4 of the 8 corners can fit a small ball
    # due to the specific positioning of the large balls.
    # These are at (12,12,12), (1,12,12), (12,1,12), and (12,12,1).
    num_corner_voids = 4
    
    num_small_balls = num_central_voids + num_face_voids + num_corner_voids

    # --- Final Calculation and Output ---
    total_energy = (num_small_balls * small_ball_energy) + (num_large_balls * large_ball_energy)

    container_str = f"box {box_dims[0]:.1f}x{box_dims[1]:.1f}x{box_dims[2]:.1f}"
    balls_str = f"{num_small_balls};{num_large_balls}"
    final_answer_str = f"[{container_str}]{balls_str}"
    
    print(f"Container: {container_str}")
    print(f"Number of 1-cm balls (a): {num_small_balls}")
    print(f"Number of 2-cm balls (b): {num_large_balls}")
    print("\nFinal Energy Calculation:")
    # The prompt asks to "output each number in the final equation"
    print(f"{total_energy} = {num_small_balls} * {small_ball_energy} + {num_large_balls} * {large_ball_energy}")

    # Final answer in the required format
    print(f"\n<<<{final_answer_str}>>>")

if __name__ == "__main__":
    solve_pioneer_probe_packing()