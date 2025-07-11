import math

def calculate_packed_spheres():
    """
    Calculates the number of spheres that can be packed into a cylinder
    using a layered heuristic approach.
    """
    # Given parameters
    r_sphere = 2.0  # Radius of the spheres
    R_cylinder = 10.0  # Radius of the cylinder
    H_cylinder = 11.0  # Height of the cylinder

    print("Problem Parameters:")
    print(f"Sphere radius (r): {r_sphere} cm")
    print(f"Cylinder radius (R): {R_cylinder} cm")
    print(f"Cylinder height (H): {H_cylinder} cm")
    print("-" * 30)

    # --- Step 1: Calculate effective dimensions for sphere centers ---
    # The space available for the centers of the spheres is smaller than the cylinder.
    h_effective = H_cylinder - 2 * r_sphere
    R_effective = R_cylinder - r_sphere

    print("Step 1: Calculate Effective Dimensions for Sphere Centers")
    print(f"Effective height for centers: H - 2r = {H_cylinder} - 2*{r_sphere} = {h_effective} cm")
    print(f"Effective radius for centers: R - r = {R_cylinder} - {r_sphere} = {R_effective} cm")
    print("-" * 30)

    # --- Step 2: Determine the number of layers ---
    # For a dense packing (like HCP), the vertical distance between the centers
    # of spheres in adjacent layers is d * sqrt(2/3), where d is the sphere diameter.
    sphere_diameter = 2 * r_sphere
    layer_height = sphere_diameter * math.sqrt(2/3)
    
    # The first layer's center is at z = r. Subsequent layers add layer_height.
    # We can fit k layers if r + (k-1)*layer_height <= H-r
    # (k-1)*layer_height <= H-2r
    # k-1 <= h_effective / layer_height
    # k <= (h_effective / layer_height) + 1
    num_layers = math.floor(h_effective / layer_height) + 1
    
    print("Step 2: Determine the Maximum Number of Layers")
    print(f"Vertical separation between layers for dense packing: 2r * sqrt(2/3) = {layer_height:.3f} cm")
    print(f"Number of layers that can fit: floor({h_effective:.1f} / {layer_height:.3f}) + 1 = {num_layers}")
    print("-" * 30)

    # --- Step 3: Determine spheres per layer ---
    # This is a 2D circle packing problem: packing circles of radius r_sphere
    # into a larger circle of radius R_effective.
    # This is equivalent to packing unit circles (r=1) into a circle of radius R_effective / r_sphere.
    ratio = R_effective / r_sphere
    
    # The optimal number of circles for a given radius ratio is a known (but hard) problem.
    # For a ratio of 4.0, the maximum number of circles is 11.
    # (From sources like packomania.com or academic papers by E. Specht).
    spheres_per_layer = 11
    
    print("Step 3: Determine Spheres per Layer")
    print(f"This is a 2D problem: Packing circles of r={r_sphere} into a circle of R={R_effective}.")
    print(f"The radius ratio is R_effective / r_sphere = {R_effective}/{r_sphere} = {ratio}")
    print(f"Based on known optimal solutions for 2D circle packing, a ratio of {ratio} can fit {spheres_per_layer} circles.")
    print("-" * 30)
    
    # --- Step 4: Calculate total spheres ---
    total_spheres = num_layers * spheres_per_layer
    
    print("Step 4: Calculate Total Number of Spheres")
    print(f"Total Spheres = (Number of Layers) * (Spheres per Layer)")
    print(f"Total Spheres = {num_layers} * {spheres_per_layer} = {total_spheres}")
    print("-" * 30)
    
    print(f"The estimated optimal number of spheres is {total_spheres}.")
    
    return total_spheres

# Run the calculation and print the final answer
final_answer = calculate_packed_spheres()
print(f"\nFinal Answer: The number of spheres is {final_answer}.")
<<<33>>>