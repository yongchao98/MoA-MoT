import math

def solve_container_design():
    """
    Calculates the minimum cost to design a container for energy balls,
    considering both box and cylinder shapes.
    """
    
    # --- Step 1: Define parameters and calculate ball cost ---
    energy_needed = 1000
    energy_per_ball = 25
    cost_per_ball = 1000
    ball_diameter = 4.0
    material_cost_per_cm2 = 200
    precision = 0.5

    num_balls = math.ceil(energy_needed / energy_per_ball)
    total_ball_cost = num_balls * cost_per_ball

    # --- Step 2: Optimize Box Container ---
    min_box_surface_area = float('inf')
    best_box_arrangement = None

    # We search for the best arrangement of balls (nx * ny * nz) that holds at least num_balls.
    # To minimize surface area, the arrangement should be as cube-like as possible.
    # We iterate through factorizations of numbers >= num_balls.
    for n_total in range(num_balls, num_balls + 25): # A reasonable search range
        for nx in range(1, int(n_total**(1/3.0)) + 2):
            if n_total % nx == 0:
                rem = n_total // nx
                for ny in range(nx, int(rem**0.5) + 2):
                    if rem % ny == 0:
                        nz = rem // ny
                        if nx * ny * nz == n_total:
                            L, W, H = nx * ball_diameter, ny * ball_diameter, nz * ball_diameter
                            area = 2 * (L*W + W*H + H*L)
                            if area < min_box_surface_area:
                                min_box_surface_area = area
                                best_box_arrangement = (nx, ny, nz)

    box_material_cost = min_box_surface_area * material_cost_per_cm2

    # --- Step 3: Optimize Cylinder Container ---
    min_cyl_surface_area = float('inf')

    # Iterate on the number of layers (nz)
    for nz in range(1, num_balls + 1):
        n_per_layer = math.ceil(num_balls / nz)
        H = nz * ball_diameter

        # Find the best 2D grid (nx_l x ny_l) for the layer to minimize the enclosing circle's radius
        min_diag_sq = float('inf')
        for nx_l in range(1, int(math.sqrt(n_per_layer)) + 2):
            # We need to pack n_per_layer balls, so we might need a grid larger than nx_l*nx_l
            ny_l = math.ceil(n_per_layer / nx_l)
            diag_sq = nx_l**2 + ny_l**2
            if diag_sq < min_diag_sq:
                min_diag_sq = diag_sq

        # Calculate container radius based on the diagonal of the layer's bounding box
        r_container = (ball_diameter / 2) * math.sqrt(min_diag_sq)
        # Round up to the nearest precision step
        R = math.ceil(r_container / precision) * precision
        
        area = 2 * math.pi * R**2 + 2 * math.pi * R * H
        if area < min_cyl_surface_area:
            min_cyl_surface_area = area
    
    cyl_material_cost = min_cyl_surface_area * material_cost_per_cm2

    # --- Step 4 & 5: Compare, Conclude, and Print Final Cost ---
    if box_material_cost < cyl_material_cost:
        chosen_design = "Box"
        chosen_area = min_box_surface_area
        chosen_material_cost = box_material_cost
    else:
        chosen_design = "Cylinder"
        chosen_area = min_cyl_surface_area
        chosen_material_cost = cyl_material_cost

    total_cost = total_ball_cost + chosen_material_cost

    print(f"Optimal design choice is a {chosen_design}.")
    print("\n--- Final Cost Calculation ---")
    print(f"Total Cost = (Number of Balls * Cost per Ball) + (Container Area * Cost per cm^2)")
    print(f"Total Cost = ({int(num_balls)} * {cost_per_ball}) + ({int(chosen_area)} * {material_cost_per_cm2})")
    print(f"Total Cost = {int(total_ball_cost)} + {int(chosen_material_cost)}")
    print(f"Total Cost = {int(total_cost)}")
    
    return total_cost

# Execute the function to get the final answer.
final_total_cost = solve_container_design()
print(f"\n<<<{int(final_total_cost)}>>>")
