import math

def solve_sphere_packing():
    """
    Calculates the number of spheres that can be packed into a cylinder
    using a layered packing heuristic.
    """
    # Given parameters
    r_s = 2.0  # Sphere radius
    R_c = 10.0 # Cylinder radius
    H_c = 11.0 # Cylinder height

    print("Problem Parameters:")
    print(f"Sphere radius (r_s): {r_s} cm")
    print(f"Cylinder radius (R_c): {R_c} cm")
    print(f"Cylinder height (H_c): {H_c} cm")
    print("-" * 30)

    # Step 1: Calculate maximum spheres in a single 2D layer.
    # This is the "circle packing in a circle" problem. The ratio of the
    # container circle radius to the inner circle radius is R_c / r_s.
    radius_ratio = R_c / r_s
    
    # We use known optimal solutions for circle packing. For a radius ratio
    # of 5, the maximum number of circles (n) is 19.
    # K_19 = 4.863... < 5, while K_20 = 5.122... > 5.
    spheres_per_layer = 19
    print("Step 1: Spheres per Layer Calculation")
    print(f"The radius ratio R_c/r_s is {R_c} / {r_s} = {radius_ratio}")
    print(f"Based on known optimal 2D packings for this ratio, the maximum number of spheres per layer is {spheres_per_layer}.")
    print("-" * 30)

    # Step 2: Calculate the number of layers that can be stacked.
    # For the densest packing (like hexagonal close-packed), the vertical
    # distance between the centers of adjacent layers is less than the sphere's diameter.
    # Vertical distance d_z = 2 * r_s * sqrt(2/3)
    vertical_separation = 2 * r_s * math.sqrt(2/3)
    
    # The total vertical space available for the centers is H_c - 2*r_s.
    available_height_for_centers = H_c - 2 * r_s
    
    # Number of layers = floor(available_height / separation) + 1
    # We can fit n separations, which means we can fit n+1 layers.
    if available_height_for_centers < 0:
        num_layers = 0
    else:
        num_layers = math.floor(available_height_for_centers / vertical_separation) + 1

    print("Step 2: Number of Layers Calculation")
    print(f"Vertical space available for sphere centers: {H_c} - 2*{r_s} = {available_height_for_centers:.2f} cm")
    print(f"Vertical separation between centers in a dense packing: {vertical_separation:.2f} cm")
    print(f"Number of layers = floor({available_height_for_centers:.2f} / {vertical_separation:.2f}) + 1 = {num_layers}")
    print("-" * 30)
    
    # Step 3: Calculate total number of spheres
    total_spheres = num_layers * spheres_per_layer

    print("Step 3: Final Calculation")
    print("The total number of spheres is the number of layers multiplied by the spheres per layer.")
    print(f"Total spheres = {num_layers} layers * {spheres_per_layer} spheres/layer = {total_spheres}")
    print("-" * 30)
    
    return total_spheres

if __name__ == '__main__':
    final_answer = solve_sphere_packing()
    print(f"\nThe estimated maximum number of spheres is: {final_answer}")
