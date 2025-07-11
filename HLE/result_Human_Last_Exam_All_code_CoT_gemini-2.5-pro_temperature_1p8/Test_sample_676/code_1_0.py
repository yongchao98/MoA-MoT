import math

# --- Problem Definition ---
BALL_RADIUS = 2.0  # cm
BALL_DIAMETER = 4.0 # cm
BALL_ENERGY = 25  # MJ
BALL_COST = 1000  # USD

MATERIAL_COST_PER_CM2 = 200  # USD

TOTAL_ENERGY_GOAL = 1000  # MJ

PRECISION = 0.5  # cm

# --- Step 1: Calculate minimum number of energy balls ---
num_balls = math.ceil(TOTAL_ENERGY_GOAL / BALL_ENERGY)
balls_total_cost = num_balls * BALL_COST

print("--- Analysis Setup ---")
print(f"Energy goal: {TOTAL_ENERGY_GOAL} MJ")
print(f"Energy per ball: {BALL_ENERGY} MJ")
print(f"Minimum number of balls required: ceil({TOTAL_ENERGY_GOAL} / {BALL_ENERGY}) = {num_balls}")
print(f"Cost of {num_balls} energy balls: {num_balls} * {BALL_COST} = ${balls_total_cost}\n")

# --- Step 2: Design and cost evaluation for a BOX container ---
print("--- Evaluating Box Container ---")
print("Finding the optimal grid packing for the balls...")

min_box_sa = float('inf')
best_box_dims = None
best_n = None

# Search for the best integer grid (nx, ny, nz) to pack num_balls
for nx in range(1, num_balls + 1):
    for ny in range(1, num_balls + 1):
        if nx * ny > num_balls:
            break
        nz = math.ceil(num_balls / (nx * ny))
        
        L = nx * BALL_DIAMETER
        W = ny * BALL_DIAMETER
        H = nz * BALL_DIAMETER
        
        sa = 2 * (L*W + L*H + W*H)
        
        if sa < min_box_sa:
            min_box_sa = sa
            best_box_dims = (L, W, H)
            best_n = (nx, ny, nz)

box_material_cost = min_box_sa * MATERIAL_COST_PER_CM2
total_box_cost = balls_total_cost + box_material_cost

print(f"Best grid found for at least {num_balls} balls: {best_n[0]} x {best_n[1]} x {best_n[2]}")
print(f"Box dimensions (L, W, H): ({best_box_dims[0]} cm, {best_box_dims[1]} cm, {best_box_dims[2]} cm)")
print(f"Box surface area: 2 * ({best_box_dims[0]}*{best_box_dims[1]} + {best_box_dims[0]}*{best_box_dims[2]} + {best_box_dims[1]}*{best_box_dims[2]}) = {min_box_sa:.2f} cm^2")
print(f"Box total cost equation: {num_balls} * {BALL_COST} + {min_box_sa:.2f} * {MATERIAL_COST_PER_CM2}")
print(f"Box total cost = ${balls_total_cost} + ${box_material_cost:.2f} = ${total_box_cost:.2f}\n")


# --- Step 3: Design and cost evaluation for a CYLINDER container ---
print("--- Evaluating Cylinder Container ---")
print("Finding the optimal layered packing for the balls...")

# R/r ratios for packing n circles in a circle, where c(n) = R_container / r_particle
# Source: http://hydra.nat.uni-magdeburg.de/packing/cci/cci.html
# c[0] is a placeholder for 1-based indexing of n circles.
packing_ratios = [0, 1.000, 2.000, 2.155, 2.414, 2.701, 3.000, 3.000, 3.305, 3.559, 3.813,
                  3.924, 3.983, 4.236, 4.328, 4.521, 4.615, 4.792, 4.864, 4.864, 5.122,
                  5.206, 5.341, 5.348, 5.567, 5.688, 5.751, 5.892, 6.000, 6.104, 6.223,
                  6.287, 6.460, 6.516, 6.641, 6.666, 6.748, 6.748, 6.756, 7.027, 7.069]

min_cyl_sa = float('inf')
best_cyl_dims = None
best_layer_config = None

for num_layers in range(1, num_balls + 1):
    balls_per_layer = math.ceil(num_balls / num_layers)
    
    container_H = num_layers * BALL_DIAMETER
    
    if balls_per_layer >= len(packing_ratios):
        continue
    
    ratio = packing_ratios[balls_per_layer]
    min_R = BALL_RADIUS * ratio
    
    # Round up radius to the nearest precision step
    container_R = math.ceil(min_R / PRECISION) * PRECISION
    
    sa = 2 * math.pi * container_R**2 + 2 * math.pi * container_R * container_H
    
    if sa < min_cyl_sa:
        min_cyl_sa = sa
        best_cyl_dims = (container_R, container_H)
        best_layer_config = (num_layers, balls_per_layer)

cyl_material_cost = min_cyl_sa * MATERIAL_COST_PER_CM2
total_cyl_cost = balls_total_cost + cyl_material_cost

(R, H) = best_cyl_dims
print(f"Best configuration: {best_layer_config[0]} layers, with up to {best_layer_config[1]} balls/layer")
print(f"Cylinder dimensions (R, H): ({R} cm, {H} cm)")
print(f"Cylinder surface area: 2*pi*{R}^2 + 2*pi*{R}*{H} = {min_cyl_sa:.2f} cm^2")
print(f"Cylinder total cost equation: {num_balls} * {BALL_COST} + {min_cyl_sa:.2f} * {MATERIAL_COST_PER_CM2}")
print(f"Cylinder total cost = ${balls_total_cost} + ${cyl_material_cost:.2f} = ${total_cyl_cost:.2f}\n")


# --- Step 4: Final comparison and result ---
print("--- Conclusion ---")
if total_box_cost < total_cyl_cost:
    print(f"The Box container is cheaper with a total cost of ${total_box_cost:.2f}.")
    final_cost = total_box_cost
    (L, W, H_box) = best_box_dims
    print(f"Final cost equation for Box: C = {num_balls} * {BALL_COST} + 2*({L}*{W} + {L}*{H_box} + {W}*{H_box}) * {MATERIAL_COST_PER_CM2}")
else:
    print(f"The Cylinder container is cheaper with a total cost of ${total_cyl_cost:.2f}.")
    final_cost = total_cyl_cost
    print(f"Final cost equation for Cylinder: C = {num_balls} * {BALL_COST} + (2*pi*{R}^2 + 2*pi*{R}*{H}) * {MATERIAL_COST_PER_CM2}")

print(f"The lowest total cost for the design is ${round(final_cost)}.")
