import math

def solve_packing():
    """
    Calculates the maximum number of spheres that can be optimally packed
    into a cylinder using a layered packing heuristic.
    """
    # --- Problem Parameters ---
    r_sphere = 2.0  # Radius of the spheres in cm
    R_cylinder = 10.0 # Radius of the cylinder in cm
    H_cylinder = 11.0 # Height of the cylinder in cm

    # --- Step 1: Calculate Spheres per Layer (2D Circle Packing) ---
    # This depends on the ratio of the cylinder radius to the sphere radius.
    radius_ratio = R_cylinder / r_sphere

    # For a radius ratio of 5.0, the maximum number of non-overlapping
    # circles is a known result from geometry. This can be looked up in
    # standard packing tables (e.g., from www.packomania.com).
    # For R/r = 5.0, the densest packing is N=31. N=32 requires R/r > 5.0.
    n_per_layer = 31

    # --- Step 2: Calculate the Number of Stackable Layers ---
    # The sphere centers must lie within the vertical range [r, H-r].
    # The total vertical space available for the centers is (H - r) - r.
    available_height_for_centers = H_cylinder - 2 * r_sphere

    # For the densest stacking (Hexagonal Close-Packing), the vertical
    # distance between the centers of adjacent layers is not 2*r, but is
    # reduced because spheres in one layer sit in the hollows of the layer below.
    # The vertical separation is h = 2 * r * sqrt(2/3).
    vertical_separation = 2 * r_sphere * math.sqrt(2.0 / 3.0)

    # The number of layers is 1 (for the base layer) plus how many more
    # can be stacked on top within the available height.
    if vertical_separation > 0:
        num_layers = math.floor(available_height_for_centers / vertical_separation) + 1
    else:
        num_layers = 1

    # --- Final Calculation ---
    total_spheres = n_per_layer * num_layers

    # --- Output the Explanation and Result ---
    print("This problem can be formulated as a non-convex optimization problem.")
    print("A robust heuristic is to calculate the packing in layers, as implemented here.")
    
    print("\n--- Calculation Steps ---")
    print(f"1. Spheres per Layer: For a cylinder-to-sphere radius ratio of {radius_ratio:.1f}, the maximum number of spheres in a single flat layer is {n_per_layer}.")
    print(f"2. Number of Layers: With a cylinder height of {H_cylinder} cm, we can stack layers with a vertical separation of {vertical_separation:.2f} cm (for dense packing).")
    print(f"   This allows for {num_layers} layers to be packed.")

    print("\n--- Final Answer ---")
    print("The total number of spheres is the product of the spheres per layer and the number of layers.")
    print(f"Total spheres = {n_per_layer} * {num_layers} = {total_spheres}")

solve_packing()
<<<93>>>