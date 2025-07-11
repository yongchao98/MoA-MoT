import math

def solve_design_problem():
    """
    Solves the container design optimization problem to find the minimum cost.
    """
    # --- Problem Parameters ---
    BALL_RADIUS = 2.0
    BALL_DIAMETER = 4.0
    BALL_ENERGY = 30  # MJ
    BALL_COST = 1000  # USD
    
    REQUIRED_ENERGY = 1000  # MJ
    MAX_SURFACE_AREA = 1000.0  # cm^2
    MATERIAL_COST_PER_CM2 = 200.0  # USD
    DIMENSION_PRECISION = 0.5  # cm

    # --- Step 1: Calculate minimum number of balls ---
    min_balls_needed = math.ceil(REQUIRED_ENERGY / BALL_ENERGY)

    best_design = {
        "cost": float('inf'),
        "shape": None,
        "surface_area": 0,
        "balls": 0,
        "dims": {}
    }

    # --- Step 2: Analyze Box Container ---
    # We assume a simple grid packing. The number of balls along each dimension
    # are nx, ny, nz. The required internal dimensions are L=nx*d, W=ny*d, H=nz*d.
    # To minimize surface area for a given number of balls, nx, ny, and nz should be
    # as close to each other as possible (cube-like).
    
    # Let's check the most cube-like arrangement for >= 34 balls: 4x3x3 = 36 balls.
    nx, ny, nz = 4, 3, 3
    num_balls_box = nx * ny * nz
    L = nx * BALL_DIAMETER  # 16 cm
    W = ny * BALL_DIAMETER  # 12 cm
    H = nz * BALL_DIAMETER  # 12 cm
    
    sa_box = 2 * (L*W + L*H + W*H)
    
    # Conclusion for the box is derived from the fact that even the most
    # efficient grid-packing configurations result in a surface area > 1000 cm^2.
    # For example, 4x3x3 balls requires a 16x12x12 cm box, with SA = 1056 cm^2.
    # Other configurations like 17x2x1 are even less efficient.
    # Therefore, we conclude no valid box design exists under these packing assumptions.

    # --- Step 3: Analyze Cylinder Container ---
    # We assume balls are stacked in layers.
    # H = num_layers * ball_diameter
    # R depends on the number of balls per layer (circle packing in a circle).
    
    # Ratios of (Container Diameter / Ball Diameter) for optimal packing of N circles.
    # Source: Wikipedia and other mathematical sources on circle packing.
    # We only need to check up to min_balls_needed, as more balls per layer
    # would lead to a very large radius and surface area.
    d_ratio_map = {
        1: 1.0, 2: 2.0, 3: 2.155, 4: 2.414, 5: 2.701, 6: 3.0, 7: 3.0,
        8: 3.304, 9: 3.613, 10: 3.813, 11: 3.923, 12: 4.029, 13: 4.236,
        14: 4.328, 15: 4.521, 16: 4.615, 17: 4.792, 18: 4.863, 19: 5.0,
        20: 5.122, 21: 5.211, 22: 5.343, 23: 5.515, 24: 5.608, 25: 5.68,
        26: 5.751, 27: 5.864, 28: 6.0, 29: 6.095, 30: 6.166, 31: 6.234,
        32: 6.32, 33: 6.425, 34: 6.478
    }

    for balls_per_layer, d_ratio in d_ratio_map.items():
        # Calculate required radius and round up to the nearest precision step
        required_radius = (d_ratio * BALL_DIAMETER) / 2.0
        container_radius = math.ceil(required_radius / DIMENSION_PRECISION) * DIMENSION_PRECISION
        
        # Calculate number of layers and height
        num_layers = math.ceil(min_balls_needed / balls_per_layer)
        container_height = num_layers * BALL_DIAMETER
        
        # Calculate total balls and surface area
        total_balls = num_layers * balls_per_layer
        surface_area = (2 * math.pi * container_radius**2) + (2 * math.pi * container_radius * container_height)
        
        if surface_area <= MAX_SURFACE_AREA:
            cost = surface_area * MATERIAL_COST_PER_CM2 + total_balls * BALL_COST
            if cost < best_design["cost"]:
                best_design["cost"] = cost
                best_design["shape"] = "Cylinder"
                best_design["surface_area"] = surface_area
                best_design["balls"] = total_balls
                best_design["dims"] = {"R": container_radius, "H": container_height}
                best_design["config"] = {"balls_per_layer": balls_per_layer, "num_layers": num_layers}

    # --- Step 4: Output the result ---
    print("--- Design Analysis ---")
    print(f"Minimum number of energy balls required: {min_balls_needed}")
    
    print("\n1. Box Container Analysis:")
    print("Assuming a simple grid packing, the most space-efficient configurations were calculated.")
    print(f"For example, to pack 36 balls (4x3x3), the required box dimensions are {L}x{W}x{H} cm.")
    print(f"This results in a surface area of {sa_box:.2f} cm^2, which exceeds the {MAX_SURFACE_AREA} cm^2 limit.")
    print("Conclusion: No valid Box design was found.")
    
    print("\n2. Cylinder Container Analysis:")
    if best_design["shape"] == "Cylinder":
        print("A valid and optimal cylinder design was found.")
        R = best_design['dims']['R']
        H = best_design['dims']['H']
        SA = best_design['surface_area']
        N_balls = best_design['balls']
        cost = best_design['cost']
        
        print("\n--- Optimal Design Details ---")
        print(f"Container Shape: {best_design['shape']}")
        print(f"Dimensions: Radius = {R:.1f} cm, Height = {H:.1f} cm")
        print(f"Configuration: {best_design['config']['balls_per_layer']} balls per layer, {best_design['config']['num_layers']} layers")
        print(f"Total Balls: {N_balls} (Energy: {N_balls * BALL_ENERGY} MJ)")
        print(f"Surface Area: {SA:.2f} cm^2")

        print("\n--- Final Cost Calculation ---")
        cost_material = SA * MATERIAL_COST_PER_CM2
        cost_balls = N_balls * BALL_COST
        
        print(f"Total Cost = (Surface Area * Cost per cm^2) + (Number of Balls * Cost per Ball)")
        print(f"Total Cost = ({SA:.2f} * {MATERIAL_COST_PER_CM2}) + ({N_balls} * {BALL_COST})")
        print(f"Total Cost = {cost_material:.2f} + {cost_balls:.2f} = {cost:.2f} USD")
        final_cost = round(cost)
    else:
        print("No valid cylinder design was found within the constraints.")
        final_cost = 0

    return final_cost


# --- Execute the function and print the final answer ---
final_answer = solve_design_problem()
print(f"\n<<<C>>>\n<<<{final_answer}>>>")
