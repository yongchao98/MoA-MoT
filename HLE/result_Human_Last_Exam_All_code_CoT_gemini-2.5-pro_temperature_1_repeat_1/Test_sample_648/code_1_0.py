import math

def solve_pioneer_mission_cost():
    """
    Solves the Pioneer mission design problem to find the minimum total cost.
    """

    # 1. Problem Parameters
    # Energy Ball properties
    ball_radius_cm = 2.0
    ball_diameter_cm = ball_radius_cm * 2
    energy_per_ball_mj = 30
    cost_per_ball_usd = 1000

    # Mission requirements
    total_energy_needed_mj = 1000

    # Container constraints
    max_surface_area_cm2 = 1000.0
    material_cost_per_cm2_usd = 200
    precision_cm = 0.5

    # 2. Calculate Ball Requirements
    num_balls_needed = math.ceil(total_energy_needed_mj / energy_per_ball_mj)
    cost_of_balls = num_balls_needed * cost_per_ball_usd

    # 3. Analyze Container Options to find the one with minimum surface area

    # The problem is to find a container that can hold `num_balls_needed` spheres
    # of diameter 4 cm, with a surface area <= 1000 cm^2.
    # We want the one with the minimum surface area.

    # Analysis of Box Container (Cubic Packing)
    # To hold `n` balls arranged in a grid of `nx x ny x nz`, the dimensions must be
    # at least L = nx*4, W = ny*4, H = nz*4.
    # To minimize surface area, we need the dimensions to be as close as possible (cube-like).
    # Let's check the most cube-like arrangement for >= 34 balls, which is 36 balls (3x3x4).
    # L=3*4=12, W=3*4=12, H=4*4=16.
    # SA_box = 2 * (12*12 + 12*16 + 12*16) = 2 * (144 + 192 + 192) = 1056 cm^2.
    # This is greater than 1000 cm^2. Less cube-like shapes will have even larger surface areas.
    # Thus, a box container with simple packing is not a feasible solution.

    # Analysis of Cylinder Container (Columnar Packing)
    # We can pack balls in vertical columns. We need to find the number of columns (`n_cols`)
    # and their arrangement to minimize the radius `R`, and the number of layers `k` to
    # determine the height `H`.
    # H = k * ball_diameter_cm
    # A highly efficient packing for circles on a plane is a hexagonal arrangement.
    # Let's test the case of 7 columns: one central column and 6 surrounding it.

    n_cols = 7  # 7 columns in a hexagonal pattern
    # The radius to contain 7 circles of radius 2cm is 6cm (center-to-outer-edge distance is 4+2=6).
    # This dimension is a multiple of 0.5cm.
    best_design = {
        "type": "Cylinder",
        "radius_cm": 6.0,
        "n_cols": n_cols
    }

    # Calculate required height
    num_layers = math.ceil(num_balls_needed / best_design["n_cols"])
    height_cm = num_layers * ball_diameter_cm
    best_design["height_cm"] = height_cm
    
    # Calculate surface area
    R = best_design["radius_cm"]
    H = best_design["height_cm"]
    surface_area = (2 * math.pi * R * H) + (2 * math.pi * R**2)
    best_design["surface_area_cm2"] = surface_area

    # Check if the design is valid
    if surface_area > max_surface_area_cm2:
        # This should not happen based on our reasoning, but it is good practice to check.
        print("No valid container design found.")
        return 0

    # 4. Calculate Final Cost
    cost_of_container = best_design["surface_area_cm2"] * material_cost_per_cm2_usd
    total_cost = cost_of_balls + cost_of_container

    # 5. Print the results and the final equation
    print("--- Optimal Design ---")
    print(f"Container Type: {best_design['type']}")
    print(f"Radius: {best_design['radius_cm']} cm")
    print(f"Height: {best_design['height_cm']} cm")
    print(f"Surface Area: {best_design['surface_area_cm2']:.2f} cm^2")
    print(f"Capacity: {best_design['n_cols'] * num_layers} balls")
    print("\n--- Cost Calculation ---")
    print(f"Energy Balls Needed: {num_balls_needed}")
    
    # Final equation and cost
    # Using the precise formula for clarity
    print("\nFinal Cost Equation:")
    print(f"Total Cost = (Balls Needed * Cost Per Ball) + (Surface Area * Material Cost)")
    
    sa_formula = f"(2 * pi * {R}^2 + 2 * pi * {R} * {H})"
    cost_formula = f"({num_balls_needed} * {cost_per_ball_usd}) + ({sa_formula} * {material_cost_per_cm2_usd})"
    
    print(f"C = {cost_formula}")
    print(f"C = {cost_of_balls} + ({best_design['surface_area_cm2']:.4f} * {material_cost_per_cm2_usd})")
    print(f"C = {cost_of_balls} + {cost_of_container:.4f}")
    print(f"C = {total_cost:.4f} USD")
    
    return total_cost

# Run the solver and print the final answer in the requested format
final_cost = solve_pioneer_mission_cost()
print(f"\n<<<C>>>\n<<<{final_cost:.2f}>>>")
