import math

def solve_container_problem():
    """
    Calculates the minimum cost to design a container for energy balls.
    """
    # Step 1: Define problem constants
    ball_radius = 2.0  # cm
    ball_diameter = ball_radius * 2.0
    ball_energy = 25  # MJ
    ball_cost = 1000  # USD
    material_cost_per_cm2 = 200  # USD
    total_energy_requirement = 1000  # MJ

    # Step 2: Calculate the minimum number of energy balls required
    min_balls_needed = total_energy_requirement / ball_energy
    print(f"Minimum number of energy balls required: {int(min_balls_needed)}")
    print("-" * 30)

    # --- Box Container Calculation ---
    print("Evaluating Box Container Design:")
    # Through optimization, the best arrangement for >= 40 balls in a grid
    # is a 2x4x5 grid.
    nl, nw, nh = 2, 4, 5
    num_balls_box = nl * nw * nh
    
    L = nl * ball_diameter
    W = nw * ball_diameter
    H = nh * ball_diameter

    surface_area_box = 2 * (L * W + L * H + W * H)
    cost_balls_box = num_balls_box * ball_cost
    cost_material_box = surface_area_box * material_cost_per_cm2
    total_cost_box = cost_balls_box + cost_material_box

    print(f"Optimal ball grid: {nl} x {nw} x {nh} = {num_balls_box} balls")
    print(f"Container dimensions: L={L:.1f} cm, W={W:.1f} cm, H={H:.1f} cm")
    print(f"Surface Area: {surface_area_box:.2f} cm^2")
    print("Total Cost (Box) = (Balls Cost) + (Material Cost)")
    print(f"Total Cost (Box) = ({num_balls_box} * ${ball_cost}) + ({surface_area_box:.2f} cm^2 * ${material_cost_per_cm2}/cm^2)")
    print(f"Total Cost (Box) = ${cost_balls_box:.2f} + ${cost_material_box:.2f} = ${total_cost_box:.2f}")
    print("-" * 30)

    # --- Cylinder Container Calculation ---
    print("Evaluating Cylinder Container Design:")
    # The optimal cylindrical packing involves stacking layers with hexagonal packing (7 balls).
    balls_per_layer_cyl = 7
    num_layers_cyl = math.ceil(min_balls_needed / balls_per_layer_cyl)
    num_balls_cyl = balls_per_layer_cyl * num_layers_cyl

    # For a hexagonal layer of 7 balls, container radius is 3 * ball_radius
    radius_cyl = 3 * ball_radius
    height_cyl = num_layers_cyl * ball_diameter

    surface_area_cyl = (2 * math.pi * radius_cyl * height_cyl) + (2 * math.pi * radius_cyl**2)
    cost_balls_cyl = num_balls_cyl * ball_cost
    cost_material_cyl = surface_area_cyl * material_cost_per_cm2
    total_cost_cyl = cost_balls_cyl + cost_material_cyl
    
    print(f"Optimal packing: {num_layers_cyl} layers of {balls_per_layer_cyl} balls = {num_balls_cyl} balls")
    print(f"Container dimensions: R={radius_cyl:.1f} cm, H={height_cyl:.1f} cm")
    print(f"Surface Area: {surface_area_cyl:.2f} cm^2")
    print("Total Cost (Cylinder) = (Balls Cost) + (Material Cost)")
    print(f"Total Cost (Cylinder) = ({num_balls_cyl} * ${ball_cost}) + ({surface_area_cyl:.2f} cm^2 * ${material_cost_per_cm2}/cm^2)")
    print(f"Total Cost (Cylinder) = ${cost_balls_cyl:.2f} + ${cost_material_cyl:.2f} = ${total_cost_cyl:.2f}")
    print("-" * 30)

    # --- Comparison and Final Answer ---
    if total_cost_cyl < total_cost_box:
        print("Conclusion: The Cylinder design is cheaper.")
        final_cost = total_cost_cyl
    else:
        print("Conclusion: The Box design is cheaper.")
        final_cost = total_cost_box
        
    print(f"\nThe minimum total cost for the project is ${final_cost:.2f}.")
    # Round to the nearest integer for the final answer format
    return round(final_cost)

if __name__ == '__main__':
    final_answer = solve_container_problem()
    print(f"<<<{final_answer}>>>")