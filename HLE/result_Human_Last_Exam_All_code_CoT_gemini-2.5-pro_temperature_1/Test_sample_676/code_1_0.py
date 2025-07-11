import math

def solve_container_problem():
    """
    Calculates the minimum cost to design a container for energy balls.
    """
    # Problem Constants
    ENERGY_REQ = 1000  # MJ
    ENERGY_PER_BALL = 25  # MJ
    COST_PER_BALL = 1000  # USD
    RADIUS_BALL = 2.0  # cm
    DIAMETER_BALL = 4.0 # cm
    COST_PER_CM2 = 200.0  # USD
    PRECISION = 0.5  # cm

    # Step 1: Calculate the number of balls and their cost
    NUM_BALLS = math.ceil(ENERGY_REQ / ENERGY_PER_BALL)
    BALLS_COST = NUM_BALLS * COST_PER_BALL

    # Step 2: Box container optimization
    # To minimize surface area for a fixed volume, dimensions should be as close as possible.
    # We need to find integer factors of 40: (nx, ny, nz) such that nx*ny*nz >= 40.
    # We will assume packing exactly 40 balls is optimal to avoid extra ball costs.
    # Factors of 40: (1,1,40), (1,2,20), (1,4,10), (1,5,8), (2,2,10), (2,4,5).
    # The set (2,4,5) is the most "cube-like".
    nx, ny, nz = 2, 4, 5
    box_L = nx * DIAMETER_BALL
    box_W = ny * DIAMETER_BALL
    box_H = nz * DIAMETER_BALL
    box_SA = 2 * (box_L * box_W + box_L * box_H + box_W * box_H)
    box_container_cost = box_SA * COST_PER_CM2
    box_total_cost = BALLS_COST + box_container_cost

    # Step 3: Cylinder container optimization
    # We need data for packing k circles in a larger circle.
    # The values are R/r ratios, where R is the containing circle's radius
    # and r is the inner circle's radius.
    pack_radii_ratio = {
        1: 1.000, 2: 2.000, 4: 2.414, 5: 2.701, 8: 3.304,
        10: 3.813, 20: 5.122, 40: 7.039
    }
    
    # Factors of 40 give pairs of (balls_per_layer, num_layers)
    factors_40 = [(1, 40), (2, 20), (4, 10), (5, 8), (8, 5), (10, 4), (20, 2), (40, 1)]
    
    min_cyl_cost = float('inf')
    best_cyl_params = {}

    for k, num_layers in factors_40:
        cyl_H = num_layers * DIAMETER_BALL
        # Calculate the minimum required radius for the container
        min_R = pack_radii_ratio[k] * RADIUS_BALL
        # Adhere to manufacturing precision by rounding up to the nearest 0.5 cm
        final_R = math.ceil(min_R / PRECISION) * PRECISION
        
        cyl_SA = 2 * math.pi * final_R**2 + 2 * math.pi * final_R * cyl_H
        container_cost = cyl_SA * COST_PER_CM2
        total_cost = BALLS_COST + container_cost
        
        if total_cost < min_cyl_cost:
            min_cyl_cost = total_cost
            best_cyl_params = {
                'cost': total_cost,
                'k': k,
                'num_layers': num_layers,
                'H': cyl_H,
                'min_R': min_R,
                'final_R': final_R,
                'SA': cyl_SA
            }
    
    # Step 4: Compare costs and print the result for the best design
    if box_total_cost < min_cyl_cost:
        final_cost = box_total_cost
        print("Optimal design: Box container")
        print(f"Number of energy balls: {NUM_BALLS}")
        print(f"Cost of Balls = {NUM_BALLS} balls * ${COST_PER_BALL:.0f}/ball = ${BALLS_COST:.0f}")
        print(f"Container configuration: {nx} x {ny} x {nz} balls")
        print(f"Container dimensions: Length = {box_L:.1f} cm, Width = {box_W:.1f} cm, Height = {box_H:.1f} cm")
        print(f"Container Surface Area = 2 * ({box_L:.1f}*{box_W:.1f} + {box_L:.1f}*{box_H:.1f} + {box_W:.1f}*{box_H:.1f}) = {box_SA:.2f} cm^2")
        print(f"Container Material Cost = {box_SA:.2f} cm^2 * ${COST_PER_CM2:.0f}/cm^2 = ${box_container_cost:.2f}")
        print(f"Total Cost = ${BALLS_COST:.2f} (Balls) + ${box_container_cost:.2f} (Container) = ${final_cost:.2f}")
    else:
        final_cost = min_cyl_cost
        params = best_cyl_params
        print("Optimal design: Cylinder container")
        print(f"Number of energy balls: {NUM_BALLS}")
        print(f"Cost of Balls = {NUM_BALLS} balls * ${COST_PER_BALL:.0f}/ball = ${BALLS_COST:.0f}")
        print(f"Container configuration: {params['num_layers']} layers of {params['k']} balls each")
        print(f"Container Height = {params['num_layers']} layers * {DIAMETER_BALL:.1f} cm/layer = {params['H']:.1f} cm")
        print(f"Container Minimum Radius = {params['min_R']:.3f} cm")
        print(f"Container Final Radius (to {PRECISION} cm precision) = {params['final_R']:.1f} cm")
        print(f"Container Surface Area = 2*pi*{params['final_R']:.1f}^2 + 2*pi*{params['final_R']:.1f}*{params['H']:.1f} = {params['SA']:.2f} cm^2")
        print(f"Container Material Cost = {params['SA']:.2f} cm^2 * ${COST_PER_CM2:.0f}/cm^2 = ${params['cost'] - BALLS_COST:.2f}")
        print(f"Total Cost = ${BALLS_COST:.2f} (Balls) + ${params['cost'] - BALLS_COST:.2f} (Container) = ${params['cost']:.2f}")

    return round(final_cost)

# Execute the function and print the final answer in the required format.
final_answer = solve_container_problem()
print(f"<<<{final_answer}>>>")
