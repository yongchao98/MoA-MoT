import math

# --- Constants ---
ENERGY_PER_BALL_MJ = 25
TOTAL_ENERGY_NEEDED_MJ = 1000
BALL_COST_USD = 1000
CONTAINER_MATERIAL_COST_USD_PER_CM2 = 200
BALL_RADIUS_CM = 2.0
BALL_DIAMETER_CM = 4.0
PRECISION_CM = 0.5

def find_best_box_design():
    """
    Calculates the minimum cost for a box-shaped container by iterating through
    possible grid arrangements of balls.
    """
    min_cost = float('inf')
    best_config = {}
    
    min_balls = math.ceil(TOTAL_ENERGY_NEEDED_MJ / ENERGY_PER_BALL_MJ)
    
    # Search for the best (nx, ny, nz) grid that holds at least min_balls.
    # To make the search efficient, we iterate nx up to min_balls,
    # and ny starting from nx to avoid duplicate shapes (e.g., 2x3x4 is same as 3x2x4).
    for nx in range(1, min_balls + 1):
        for ny in range(nx, min_balls + 1):
            # Calculate the number of layers (nz) needed
            if nx * ny == 0: continue
            nz = math.ceil(min_balls / (nx * ny))
            
            num_balls = nx * ny * nz
            
            # Dimensions of the box
            length = nx * BALL_DIAMETER_CM
            width = ny * BALL_DIAMETER_CM
            height = nz * BALL_DIAMETER_CM
            
            surface_area = 2 * (length * width + length * height + width * height)
            
            balls_cost = num_balls * BALL_COST_USD
            container_cost = surface_area * CONTAINER_MATERIAL_COST_USD_PER_CM2
            total_cost = balls_cost + container_cost
            
            if total_cost < min_cost:
                min_cost = total_cost
                best_config = {
                    "type": "Box", "cost": total_cost, "balls_cost": balls_cost,
                    "container_cost": container_cost, "num_balls": num_balls,
                    "nx": nx, "ny": ny, "nz": nz, "L": length, "W": width, 
                    "H": height, "SA": surface_area
                }
            # Optimization: if nx*ny is large, nz will be 1, and increasing ny further
            # will only make the shape less cube-like and more expensive.
            if nx * ny > min_balls:
                break
    return best_config

def find_best_cylinder_design():
    """
    Calculates the minimum cost for a cylinder-shaped container by testing
    different numbers of balls per layer.
    """
    min_cost = float('inf')
    best_config = {}

    min_balls = math.ceil(TOTAL_ENERGY_NEEDED_MJ / ENERGY_PER_BALL_MJ)
    
    # Data for packing k circles of radius r in a larger circle.
    # Value is the ratio R_min/r for the best known packing.
    packing_ratios = {
        1: 1.0, 2: 2.0, 3: 2.155, 4: 2.414, 5: 2.701, 6: 3.0, 7: 3.0, 
        8: 3.305, 9: 3.532, 10: 3.813, 11: 3.924, 12: 4.029, 13: 4.236, 
        14: 4.328, 15: 4.521, 16: 4.615, 17: 4.792, 18: 4.864, 19: 4.864
    }

    # Iterate through k (number of balls per layer)
    for k in range(1, min_balls + 1):
        if k not in packing_ratios:
            # Stop if we run out of packing data, assuming best is for small k
            break

        min_radius = packing_ratios[k] * BALL_RADIUS_CM
        container_radius = math.ceil(min_radius / PRECISION_CM) * PRECISION_CM
        num_layers = math.ceil(min_balls / k)
        container_height = num_layers * BALL_DIAMETER_CM
        num_balls = k * num_layers

        surface_area = (2 * math.pi * container_radius**2) + (2 * math.pi * container_radius * container_height)
        
        balls_cost = num_balls * BALL_COST_USD
        container_cost = surface_area * CONTAINER_MATERIAL_COST_USD_PER_CM2
        total_cost = balls_cost + container_cost
        
        if total_cost < min_cost:
            min_cost = total_cost
            best_config = {
                "type": "Cylinder", "cost": total_cost, "balls_cost": balls_cost,
                "container_cost": container_cost, "num_balls": num_balls, "k": k, 
                "n_l": num_layers, "R": container_radius, "H": container_height, "SA": surface_area
            }
    return best_config
    
# --- Main Execution ---
box_design = find_best_box_design()
cylinder_design = find_best_cylinder_design()

if box_design.get("cost", float('inf')) < cylinder_design.get("cost", float('inf')):
    best_design = box_design
else:
    best_design = cylinder_design

print("Comparing the two possible designs, the optimal design is chosen to minimize total cost.\n")
print(f"Final Recommended Design: {best_design['type']}")
print("---------------------------------------------")
print("The final cost C is calculated based on this design.")
print("Here is the breakdown of the final cost:\n")

# Output the detailed equation for the final cost
print(f"1. Number of Energy Balls (N): {best_design['num_balls']}")
print(f"   Cost of Balls = {best_design['num_balls']} balls * ${BALL_COST_USD}/ball = ${best_design['balls_cost']:.2f}\n")

if best_design['type'] == 'Cylinder':
    R, H, SA = best_design['R'], best_design['H'], best_design['SA']
    print("2. Container Details (Cylinder):")
    print(f"   - Optimal balls per layer: {best_design['k']}")
    print(f"   - Number of layers: {best_design['n_l']}")
    print(f"   - Container Radius (R): {R:.1f} cm")
    print(f"   - Container Height (H): {H:.1f} cm")
    print(f"   - Surface Area (SA) = 2 * pi * R * (R + H)")
    print(f"   - SA = 2 * {math.pi:.4f} * {R:.1f} * ({R:.1f} + {H:.1f}) = {SA:.2f} cm^2")
else:
    L, W, H, SA = best_design['L'], best_design['W'], best_design['H'], best_design['SA']
    print("2. Container Details (Box):")
    print(f"   - Optimal grid of balls: {best_design['nx']} x {best_design['ny']} x {best_design['nz']}")
    print(f"   - Container Length (L): {L:.1f} cm")
    print(f"   - Container Width (W): {W:.1f} cm")
    print(f"   - Container Height (H): {H:.1f} cm")
    print(f"   - Surface Area (SA) = 2 * (L*W + L*H + W*H)")
    print(f"   - SA = 2 * ({L:.1f}*{W:.1f} + {L:.1f}*{H:.1f} + {W:.1f}*{H:.1f}) = {SA:.2f} cm^2")

print(f"\n   Cost of Container = {SA:.2f} cm^2 * ${CONTAINER_MATERIAL_COST_USD_PER_CM2}/cm^2 = ${best_design['container_cost']:.2f}\n")

print("3. Total Cost Calculation:")
print(f"   Total Cost (C) = Cost of Balls + Cost of Container")
print(f"   C = ${best_design['balls_cost']:.2f} + ${best_design['container_cost']:.2f}")
print(f"   C = ${best_design['cost']:.2f}")

final_cost_value = best_design.get('cost', 0)
# Final answer wrapper
# print(f"\n<<<C = {final_cost_value}>>>") # Per instruction, final answer is outside code block