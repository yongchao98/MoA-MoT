import math

def solve_pioneer_probe_design():
    """
    Solves the Pioneer probe container design problem by finding the
    minimum cost cylindrical container that meets all specifications.
    """
    
    # Cost and physical parameters
    MIN_REQUIRED_ENERGY = 1000  # MJ
    ENERGY_PER_BALL = 30  # MJ
    BALL_RADIUS = 2  # cm
    BALL_DIAMETER = 4 # cm
    BALL_COST = 1000  # usd
    MATERIAL_COST_PER_CM2 = 200  # usd/cm^2
    MAX_SURFACE_AREA = 1000  # cm^2
    DIMENSION_PRECISION = 0.5 # cm
    RADIUS_PRECISION = DIMENSION_PRECISION / 2.0

    # Minimum number of balls required
    min_balls = math.ceil(MIN_REQUIRED_ENERGY / ENERGY_PER_BALL)

    # We pre-calculate that a box is not feasible. The most compact
    # box holding 36 balls (4x3x3 packing) has SA=1056 cm^2, > 1000 cm^2.
    # We focus on the cylindrical container.
    
    # Data for packing n unit circles into a larger circle.
    # The value is the radius of the container circle for unit-radius inner circles.
    # Source: http://hydra.nat.uni-magdeburg.de/packing/cpc/
    # We multiply by BALL_RADIUS (2cm) to get the required container radius.
    k_n = {
        1: 1.0, 2: 2.0, 3: 2.1547, 4: 2.4142, 5: 2.7013, 6: 3.0, 7: 3.0,
        8: 3.3047, 9: 3.5023, 10: 3.8130, 11: 3.9238, 12: 3.9829, 13: 4.0849,
        14: 4.2105, 15: 4.3283, 16: 4.4230, 17: 4.5, 18: 4.6154, 19: 4.6930,
        20: 4.7924, 21: 4.8637, 22: 4.9392, 23: 5.0118, 24: 5.0772, 25: 5.1613,
        26: 5.2349, 27: 5.3094, 28: 5.3615, 29: 5.4211, 30: 5.5, 31: 5.5562,
        32: 5.6174, 33: 5.6888, 34: 5.7533
    }
    
    min_total_cost = float('inf')
    best_design = None

    # Iterate through possible number of balls per layer (n_L)
    for n_L in range(1, int(min_balls) + 1):
        # Calculate required number of layers
        n_H = math.ceil(min_balls / n_L)
        
        # Total number of balls for this configuration
        num_balls = n_L * n_H
        
        # Calculate container dimensions
        container_H = n_H * BALL_DIAMETER
        
        # Calculate required container radius
        min_R = k_n[n_L] * BALL_RADIUS
        
        # Adjust radius for manufacturing precision (must be multiple of 0.25 cm)
        container_R = math.ceil(min_R / RADIUS_PRECISION) * RADIUS_PRECISION

        # Calculate surface area
        surface_area = 2 * math.pi * container_R * (container_R + container_H)
        
        if surface_area <= MAX_SURFACE_AREA:
            # Calculate total cost
            cost_balls = num_balls * BALL_COST
            cost_container = surface_area * MATERIAL_COST_PER_CM2
            total_cost = cost_balls + cost_container

            if total_cost < min_total_cost:
                min_total_cost = total_cost
                best_design = {
                    "cost": total_cost,
                    "balls": num_balls,
                    "sa": surface_area,
                    "n_L": n_L,
                    "n_H": n_H,
                    "R": container_R,
                    "H": container_H
                }

    if best_design:
        cost_balls = best_design["balls"] * BALL_COST
        cost_container = best_design["sa"] * MATERIAL_COST_PER_CM2
        print(f"Optimal design found: Cylindrical Container")
        print(f"Configuration: {best_design['n_L']} balls/layer, {best_design['n_H']} layers")
        print(f"Dimensions: Radius = {best_design['R']:.2f} cm, Height = {best_design['H']:.2f} cm")
        print(f"Total Balls: {best_design['balls']}, Surface Area: {best_design['sa']:.2f} cm^2")
        print(f"Final Cost Equation:")
        print(f"Total Cost = ({best_design['balls']} balls * {BALL_COST} usd/ball) + ({best_design['sa']:.2f} cm^2 * {MATERIAL_COST_PER_CM2} usd/cm^2)")
        print(f"Total Cost = {cost_balls} + {cost_container:.2f} = {best_design['cost']:.2f} usd")
        print(f"<<<{best_design['cost']:.2f}>>>")
    else:
        print("No feasible design found that meets all constraints.")
        print("<<<0>>>")

solve_pioneer_probe_design()