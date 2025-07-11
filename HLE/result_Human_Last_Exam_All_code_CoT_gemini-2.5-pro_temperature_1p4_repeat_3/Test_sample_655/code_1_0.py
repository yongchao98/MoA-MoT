import math

def solve_pioneer_probe_packing():
    """
    This function calculates the optimal packing of energy balls for the Pioneer probe
    based on a structured packing strategy in the best-found container shape.
    """

    # --- Configuration based on analysis ---
    # Container: A 13.0cm x 13.0cm x 13.0cm box.
    # SA = 6 * 13^2 = 1014 cm^2 (which is <= 1050 cm^2).
    container_dims = {'L': 13.0, 'W': 13.0, 'H': 13.0}
    container_half_dims = {k: v / 2.0 for k, v in container_dims.items()}

    # Strategy:
    # 1. Pack a 3x3x3 grid of 2-cm radius balls.
    # 2. Pack 1-cm radius balls in the interstitial sites.

    big_balls_centers = []
    big_ball_radius = 2.0

    # Place the 2-cm radius balls on a 4cm grid centered at the origin
    for i in [-1, 0, 1]:
        for j in [-1, 0, 1]:
            for k in [-1, 0, 1]:
                center_pos = (4.0 * i, 4.0 * j, 4.0 * k)
                big_balls_centers.append(center_pos)

    # Verify all big balls are inside the container
    for center in big_balls_centers:
        # Farthest point of a ball from origin is its center + radius
        if not (abs(center[0]) + big_ball_radius <= container_half_dims['L'] and
                abs(center[1]) + big_ball_radius <= container_half_dims['W'] and
                abs(center[2]) + big_ball_radius <= container_half_dims['H']):
            # This should not happen with our chosen dimensions
            raise ValueError(f"Big ball at {center} is outside the container.")
    
    # After placing big balls, create a list of all packed balls so far
    packed_balls = [(c, big_ball_radius) for c in big_balls_centers]

    small_balls_centers = []
    small_ball_radius = 1.0

    # Try to place 1-cm radius balls in the center of the cubelets
    # formed by the big balls' centers.
    for i in [-1, 1]:
        for j in [-1, 1]:
            for k in [-1, 1]:
                center_pos = (2.0 * i, 2.0 * j, 2.0 * k)
                
                # Check if this ball fits inside the container
                is_inside = (abs(center_pos[0]) + small_ball_radius <= container_half_dims['L'] and
                             abs(center_pos[1]) + small_ball_radius <= container_half_dims['W'] and
                             abs(center_pos[2]) + small_ball_radius <= container_half_dims['H'])
                
                if not is_inside:
                    continue
                
                # Check for collision with any already packed ball
                can_place = True
                for p_center, p_radius in packed_balls:
                    dist_sq = sum((c1 - c2)**2 for c1, c2 in zip(center_pos, p_center))
                    min_dist_sq = (small_ball_radius + p_radius)**2
                    if dist_sq < min_dist_sq - 1e-9:  # Using tolerance for float comparison
                        can_place = False
                        break
                
                if can_place:
                    small_balls_centers.append(center_pos)

    num_small_balls = len(small_balls_centers)
    num_big_balls = len(big_balls_centers)
    
    # Format the final answer string
    container_desc = f"box {container_dims['L']}x{container_dims['W']}x{container_dims['H']}"
    final_answer_string = f"[{container_desc}]{num_small_balls};{num_big_balls}"
    
    # Print the equation representing the final answer components
    print(f"Container description: {container_desc}")
    print(f"Number of 1-cm balls (a): {num_small_balls}")
    print(f"Number of 2-cm balls (b): {num_big_balls}")
    print(f"Final answer format: [C]a;b")
    print(f"Result: {final_answer_string}")

# Execute the function
solve_pioneer_probe_packing()

# The final answer in the required format
final_answer = "[box 13.0x13.0x13.0]8;27"
print(f"\n<<<__{final_answer}__>>>")

# Silently correct the output format to not have double underscores.
final_answer = final_answer.replace("__", "")
# print(f"\n<<<{final_answer}>>>")