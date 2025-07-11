import math

def solve_design_problem():
    """
    Finds the optimal container design to minimize total cost.
    """
    # Problem constants
    ENERGY_REQUIRED = 1000  # MJ
    ENERGY_PER_BALL = 30  # MJ
    COST_PER_BALL = 1000  # usd
    MAX_SURFACE_AREA = 1000  # cm^2
    COST_PER_CM2 = 200  # usd
    BALL_DIAMETER = 4.0  # cm
    DIMENSION_PRECISION = 0.5  # cm

    # Step 1: Calculate minimum number of balls
    min_balls = math.ceil(ENERGY_REQUIRED / ENERGY_PER_BALL)

    # Step 2: Calculate the required packing volume for these balls
    # We assume each ball occupies a cubic space of side = ball_diameter
    packing_volume_per_ball = BALL_DIAMETER ** 3
    required_volume = min_balls * packing_volume_per_ball

    min_total_cost = float('inf')
    best_design = None

    # Step 3: Search for the optimal cylinder dimensions
    # A box is not feasible as its max volume for SA<=1000 is too small.
    # Search range for radius R (as a multiple of 0.5cm)
    # R must be at least 2.0. Max R is when H=0, 2*pi*R^2=1000 -> R~12.6
    # So r_steps goes from 4 (for R=2.0) to 25 (for R=12.5)
    for r_steps in range(int(BALL_DIAMETER / 2 / DIMENSION_PRECISION), 26):
        R = r_steps * DIMENSION_PRECISION
        
        # Search range for height H. Max H is when R is smallest (R=2), H ~ 77.5
        # h_steps goes up to 160 (for H=80.0)
        for h_steps in range(int(BALL_DIAMETER / DIMENSION_PRECISION), 161):
            H = h_steps * DIMENSION_PRECISION

            surface_area = 2 * math.pi * R * (R + H)
            
            # Check surface area constraint
            if surface_area <= MAX_SURFACE_AREA:
                volume = math.pi * (R ** 2) * H
                
                # Check if it can hold the minimum number of balls
                if volume >= required_volume:
                    
                    # This design is valid for the minimum number of balls (34)
                    num_balls_to_use = min_balls
                    
                    container_cost = surface_area * COST_PER_CM2
                    balls_cost = num_balls_to_use * COST_PER_BALL
                    total_cost = container_cost + balls_cost

                    if total_cost < min_total_cost:
                        min_total_cost = total_cost
                        best_design = {
                            "type": "Cylinder",
                            "R": R,
                            "H": H,
                            "SA": surface_area,
                            "balls": num_balls_to_use,
                            "cost": total_cost
                        }

    if best_design:
        # Step 4: Output the result
        R_opt = best_design['R']
        H_opt = best_design['H']
        SA_opt = best_design['SA']
        N_opt = best_design['balls']
        C_opt = best_design['cost']
        
        print("Optimal Design Found:")
        print(f"Container Type: {best_design['type']}")
        print(f"Radius: {R_opt} cm")
        print(f"Height: {H_opt} cm")
        print(f"Number of Energy Balls: {N_opt}")
        print(f"Surface Area: {SA_opt:.2f} cm^2")
        print(f"Total Cost: ${C_opt:.2f}")
        print("\nFinal Cost Calculation:")
        # The prompt requires showing the final equation with numbers
        print(f"{C_opt:.2f} = ({SA_opt:.4f} * {COST_PER_CM2}) + ({N_opt} * {COST_PER_BALL})")

        return C_opt
    else:
        print("No feasible solution found.")
        return 0

# Run the solver and get the final cost C
C = solve_design_problem()

# The final answer format as requested by the prompt
print(f"\n<<<C={C:.2f}>>>")