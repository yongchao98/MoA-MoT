import math

def solve_design_problem():
    """
    Calculates the minimum cost for a container and energy balls
    to meet the mission requirements.
    """

    # --- Problem Parameters ---
    total_energy_req_mj = 1000
    energy_per_ball_mj = 30
    cost_per_ball_usd = 1000
    
    ball_radius_cm = 2
    ball_diameter_cm = ball_radius_cm * 2
    
    max_surface_area_cm2 = 1000
    material_cost_per_cm2 = 200
    
    precision_cm = 0.5

    # Step 1: Calculate minimum number of balls needed
    min_balls_needed = math.ceil(total_energy_req_mj / energy_per_ball_mj)
    
    print(f"--- Step 1: Analyze Energy Requirement ---")
    print(f"Total energy required: {total_energy_req_mj} MJ")
    print(f"Energy per ball: {energy_per_ball_mj} MJ")
    print(f"Minimum number of balls = ceil({total_energy_req_mj} / {energy_per_ball_mj}) = {min_balls_needed}\n")

    # Step 2: Find the optimal container configuration
    # Through analysis, simple cubic packing is not feasible due to the surface area constraint.
    # We must assume a denser packing method. We will use Face-Centered Cubic (FCC) packing.
    # The most cube-like configuration is best for minimizing surface area.
    # We will test a 3x3x4 arrangement, which can hold up to 36 balls.
    
    n_L, n_W, n_H = 3, 3, 4
    num_balls_in_container = min_balls_needed

    print(f"--- Step 2: Design the Container (Box) ---")
    print(f"Assuming an efficient FCC (Face-Centered Cubic) packing arrangement.")
    print(f"The most balanced configuration to hold at least {min_balls_needed} balls is a {n_L}x{n_W}x{n_H} grid.\n")

    # Step 3: Calculate container dimensions
    # For a box with FCC packing, the base dimensions are standard.
    # The height is compressed compared to simple cubic packing.
    # Height formula for n_H layers: H = (n_H - 1) * d * sqrt(2/3) + d
    
    l_req = n_L * ball_diameter_cm
    w_req = n_W * ball_diameter_cm
    # The height for FCC packing of n_H layers is smaller than n_H * ball_diameter
    h_req = (n_H - 1) * ball_diameter_cm * math.sqrt(2/3) + ball_diameter_cm

    # Apply manufacturing precision
    L = math.ceil(l_req / precision_cm) * precision_cm
    W = math.ceil(w_req / precision_cm) * precision_cm
    H = math.ceil(h_req / precision_cm) * precision_cm

    print(f"--- Step 3: Calculate Final Dimensions ---")
    print(f"Required internal dimensions for a {n_L}x{n_W}x{n_H} FCC grid:")
    print(f"Length = {n_L} * {ball_diameter_cm} = {l_req:.3f} cm")
    print(f"Width = {n_W} * {ball_diameter_cm} = {w_req:.3f} cm")
    print(f"Height = ({n_H}-1)*{ball_diameter_cm}*sqrt(2/3) + {ball_diameter_cm} = {h_req:.3f} cm")
    print(f"\nDimensions after applying {precision_cm} cm manufacturing precision:")
    print(f"Final Length = {L:.1f} cm")
    print(f"Final Width = {W:.1f} cm")
    print(f"Final Height = {H:.1f} cm\n")

    # Step 4: Calculate Surface Area and check constraint
    surface_area = 2 * (L * W + L * H + W * H)

    print(f"--- Step 4: Calculate Surface Area ---")
    print(f"Surface Area = 2 * (L*W + L*H + W*H)")
    print(f"Surface Area = 2 * ({L:.1f}*{W:.1f} + {L:.1f}*{H:.1f} + {W:.1f}*{H:.1f})")
    print(f"Surface Area = 2 * ({L*W:.2f} + {L*H:.2f} + {W*H:.2f}) = {surface_area:.2f} cm^2")
    
    if surface_area > max_surface_area_cm2:
        print(f"Error: Calculated surface area {surface_area:.2f} exceeds the maximum of {max_surface_area_cm2}.")
        return 0
    else:
        print(f"The surface area is within the {max_surface_area_cm2} cm^2 limit.\n")

    # Step 5: Calculate total cost
    cost_material = surface_area * material_cost_per_cm2
    cost_balls = num_balls_in_container * cost_per_ball_usd
    total_cost = cost_material + cost_balls
    
    print(f"--- Step 5: Calculate Total Cost ---")
    print("Total Cost = (Surface Area * Material Cost) + (Number of Balls * Ball Cost)")
    print(f"Material Cost = {surface_area:.2f} cm^2 * ${material_cost_per_cm2}/cm^2 = ${cost_material:,.0f}")
    print(f"Energy Ball Cost = {num_balls_in_container} balls * ${cost_per_ball_usd}/ball = ${cost_balls:,.0f}")
    print(f"Total Cost = ${cost_material:,.0f} + ${cost_balls:,.0f} = ${total_cost:,.0f}")

    return total_cost

# --- Execute the solution ---
final_cost = solve_design_problem()
# The final answer format
if final_cost is not None:
    print(f"<<<{final_cost:.0f}>>>")
