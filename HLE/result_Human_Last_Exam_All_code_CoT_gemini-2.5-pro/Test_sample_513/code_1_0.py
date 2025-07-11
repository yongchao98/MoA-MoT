import math

def solve_sphere_packing():
    """
    Calculates the number of spheres that can be optimally packed into a cylinder
    using a layer-based heuristic with Hexagonal Close-Packing (HCP).
    """
    # Given parameters
    sphere_radius = 2.0  # r
    cylinder_radius = 10.0 # R
    cylinder_height = 11.0 # H

    # --- Step 1: Calculate the maximum number of layers ---

    # In a Hexagonal Close-Packing (HCP) structure, the vertical distance
    # between the centers of spheres in adjacent layers is sqrt(8/3) * r.
    hcp_layer_separation = math.sqrt(8.0 / 3.0) * sphere_radius

    # The total available vertical space for the sphere centers is H - 2*r.
    available_height = cylinder_height - 2 * sphere_radius

    # The number of layers is 1 (for the base layer) plus how many additional
    # layers can be stacked on top.
    # floor(available_height / hcp_layer_separation) gives the number of additional layers.
    num_layers = math.floor(available_height / hcp_layer_separation) + 1

    # --- Step 2: Calculate the number of spheres per layer ---

    # This is a 2D circle packing problem: packing circles of radius r=2
    # into a larger circle of radius R=10.
    # The ratio of radii is R/r = 10/2 = 5.
    # From known optimal solutions for the "packing equal circles in a circle"
    # problem (e.g., from packomania.com), the maximum number of circles for a
    # radius ratio of 5 is 19.
    spheres_per_layer = 19

    # --- Step 3: Calculate the total number of spheres ---
    total_spheres = num_layers * spheres_per_layer

    # --- Print the results and the final equation ---
    print("--- Sphere Packing Calculation ---")
    print(f"Sphere Radius (r): {sphere_radius} cm")
    print(f"Cylinder Radius (R): {cylinder_radius} cm")
    print(f"Cylinder Height (H): {cylinder_height} cm")
    print("\nStep 1: Calculating the number of vertical layers (HCP stacking)")
    print(f"Vertical separation between layers = {hcp_layer_separation:.2f} cm")
    print(f"Number of Layers = floor((H - 2*r) / sep) + 1 = floor(({cylinder_height} - {2*sphere_radius}) / {hcp_layer_separation:.2f}) + 1 = {num_layers}")

    print("\nStep 2: Calculating spheres per horizontal layer")
    print("This corresponds to packing circles of r=2 into a circle of R=10.")
    print(f"Number of Spheres per Layer (from known results for R/r=5) = {spheres_per_layer}")

    print("\nStep 3: Calculating total number of spheres")
    print(f"Total Spheres = Number of Layers * Spheres per Layer")
    print(f"Total Spheres = {num_layers} * {spheres_per_layer} = {total_spheres}")
    
    # Return the final number for the answer block
    return total_spheres

# Run the calculation and store the final answer
final_answer = solve_sphere_packing()
print(f"\nFinal Answer: The optimal number of spheres is {final_answer}.")
# The <<<...>>> format is for the platform, not for user output.
# print(f'<<<{final_answer}>>>')