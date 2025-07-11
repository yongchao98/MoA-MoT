import math

def solve():
    """
    This function analyzes the energy packing problem and prints the optimal solution found.
    """

    # --- Problem Parameters ---
    MAX_SURFACE_AREA = 1050.0  # cm^2
    
    # Energy ball specifications
    # Ball 1 (small): 1cm radius, 1 MJ energy
    r1 = 1.0
    d1 = 2.0
    e1 = 1
    # Ball 2 (large): 2cm radius, 10 MJ energy
    r2 = 2.0
    d2 = 4.0
    e2 = 10

    # --- Container Candidate: A specific Rectangular Box ---
    l, w, h = 12.0, 12.0, 15.5
    container_description = f"box {l}x{w}x{h}"

    # --- Verification Step 1: Surface Area Constraint ---
    surface_area = 2 * (l*w + l*h + w*h)
    print("--- Analysis of the Proposed Box Container Solution ---")
    print(f"Container dimensions: {l}cm x {w}cm x {h}cm")
    print(f"Calculated surface area: 2 * ({l}*{w} + {l}*{h} + {w}*{h}) = {surface_area:.2f} cm^2")
    if surface_area > MAX_SURFACE_AREA:
        print(f"Error: Surface area {surface_area} exceeds the maximum of {MAX_SURFACE_AREA}")
        return
    print(f"Surface area is within the allowed maximum of {MAX_SURFACE_AREA} cm^2.\n")


    # --- Verification Step 2: Packing Calculation ---
    # First, pack the high-energy (large) balls.
    # A 3x3x3 grid of 2-cm radius balls fits perfectly within a 12x12x12 cm space.
    nx_b = l / d2
    ny_b = w / d2
    nz_b_in_cube = l / d2 
    
    num_big_balls = math.floor(nx_b) * math.floor(ny_b) * math.floor(nz_b_in_cube)
    
    print("--- Packing Calculation ---")
    print(f"1. Large Balls (radius {r2} cm, energy {e2} MJ):")
    print(f"A simple cubic lattice of {math.floor(nx_b)}x{math.floor(ny_b)}x{math.floor(nz_b_in_cube)} is placed in the {l}x{w}x{l} cm part of the box.")
    print(f"Number of 2-cm radius balls (b): {num_big_balls}")
    
    # Second, pack low-energy (small) balls in the remaining space.
    # This space consists of two types of voids:
    # 1. Interstitial voids within the 3x3x3 lattice of large balls.
    # 2. An empty slab on top measuring 12x12x(15.5-12) = 12x12x3.5 cm.
    
    print(f"\n2. Small Balls (radius {r1} cm, energy {e1} MJ):")
    # 1. Packing in interstitial voids
    # A 3x3x3 lattice has (3-1)x(3-1)x(3-1) = 8 central voids large enough for a 1-cm ball.
    num_small_balls_voids = (math.floor(nx_b) - 1) * (math.floor(ny_b) - 1) * (math.floor(nz_b_in_cube) - 1)
    print(f" - {num_small_balls_voids} balls are placed in the interstitial voids of the large ball lattice.")

    # 2. Packing in the empty slab
    slab_l, slab_w, slab_h = l, w, h - l
    nx_s = slab_l / d1
    ny_s = slab_w / d1
    nz_s = slab_h / d1
    num_small_balls_slab = math.floor(nx_s) * math.floor(ny_s) * math.floor(nz_s)
    print(f" - {num_small_balls_slab} balls are packed in a {math.floor(nx_s)}x{math.floor(ny_s)}x{math.floor(nz_s)} grid in the remaining slab ({slab_l}x{slab_w}x{slab_h} cm).")

    num_small_balls = num_small_balls_voids + num_small_balls_slab
    print(f"Total number of 1-cm radius balls (a): {num_small_balls}")
    
    # --- Final Calculation & Result ---
    total_energy = num_small_balls * e1 + num_big_balls * e2
    print("\n--- Summary ---")
    print(f"Container: {container_description}")
    print(f"Number of 1-cm balls (a): {num_small_balls}")
    print(f"Number of 2-cm balls (b): {num_big_balls}")
    print(f"Total potential energy equation: {num_small_balls} * {e1} MJ + {num_big_balls} * {e2} MJ = {total_energy} MJ")

    final_answer_format = f"[{container_description}]{num_small_balls};{num_big_balls}"
    
    # The final answer is printed in the required format.
    print(f"<<<{final_answer_format}>>>")

if __name__ == '__main__':
    solve()