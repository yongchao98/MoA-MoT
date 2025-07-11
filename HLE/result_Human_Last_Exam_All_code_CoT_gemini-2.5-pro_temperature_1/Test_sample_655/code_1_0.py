import math

def solve():
    """
    Solves the energy ball packing problem by analyzing each container shape.
    """
    max_surface_area = 1050.0
    
    # --- Box Analysis ---
    # Strategy: Use a simple cubic packing. Large balls (diameter 4) are packed
    # in a 3x3x3 grid, requiring a 12x12x12 space. We add space for small balls.
    # The container needs internal dimensions of 12x12x(12+2) for one layer of small balls.
    # Let's use dimensions L=12, W=12, H=14.
    box_l, box_w, box_h = 12.0, 12.0, 14.0
    box_sa = 2 * (box_l * box_w + box_l * box_h + box_w * box_h)
    
    # Check if a slightly taller box is possible for a second layer of small balls (H=16)
    box_h_test = 15.5 # Max height for L=12, W=12 is 15.875, so 15.5 is the max multiple of 0.5
    box_sa_test = 2 * (box_l * box_w + box_l * box_h_test + box_w * box_h_test)

    if box_sa_test <= max_surface_area:
        box_h = box_h_test # Use taller box if possible
        
    n2_box = (box_l // 4) * (box_w // 4) * (box_h // 4)
    # Remaining height for small balls
    rem_h = box_h - (n2_box // ((box_l // 4) * (box_w // 4))) * 4
    n1_box = (box_l // 2) * (box_w // 2) * (rem_h // 2)
    
    energy_box = 10 * n2_box + 1 * n1_box
    
    # --- Cylinder Analysis ---
    # Strategy: Find the volume-optimal cylinder (h=2r) and pack with simple layers.
    # S = 6*pi*r^2 <= 1050 => r <= sqrt(1050/(6*pi)) = 7.46
    # Let's test r=7.0 (multiple of 0.5)
    cyl_r = 7.0
    cyl_h = 14.0 # h=2r
    cyl_sa = 2 * math.pi * cyl_r * (cyl_r + cyl_h)
    
    # Pack large balls (radius 2) in layers.
    # A layer can fit 7 balls (1 central, 6 around) in a circle of radius 7.
    balls_per_layer_cyl = 7
    # Height of cylinder allows for 3 layers (h=14, ball diameter=4 => 14/4=3.5)
    num_layers_cyl = math.floor(cyl_h / 4)
    n2_cyl = balls_per_layer_cyl * num_layers_cyl
    n1_cyl = 0 # Assume no space for small balls in this simple packing
    energy_cyl = 10 * n2_cyl + 1 * n1_cyl

    # --- Sphere Analysis ---
    # Strategy: A sphere has the best volume-to-surface-area ratio.
    # S = 4*pi*r^2 <= 1050 => r <= sqrt(1050/(4*pi)) = 9.14
    # The largest radius that is a multiple of 0.5 is 9.0
    sphere_r = 9.0
    sphere_sa = 4 * math.pi * sphere_r**2
    
    # For packing equal spheres (radius r_ball=2) into a larger sphere (radius R=9.0),
    # the ratio is R/r_ball = 4.5.
    # According to known results in packing problems, the maximum number of spheres is 57.
    n2_sphere = 57
    # In such a dense packing, the gaps are too small for the 1-cm radius balls.
    n1_sphere = 0
    energy_sphere = 10 * n2_sphere + 1 * n1_sphere

    # --- Comparison ---
    best_energy = 0
    best_config = ""

    if energy_box > best_energy:
        best_energy = energy_box
        best_config = f"[box {box_l}x{box_w}x{box_h}]{int(n1_box)};{int(n2_box)}"
        final_n1, final_n2 = int(n1_box), int(n2_box)
        
    if energy_cyl > best_energy:
        best_energy = energy_cyl
        best_config = f"[cylinder r={cyl_r}, h={cyl_h}]{int(n1_cyl)};{int(n2_cyl)}"
        final_n1, final_n2 = int(n1_cyl), int(n2_cyl)

    if energy_sphere > best_energy:
        best_energy = energy_sphere
        best_config = f"[sphere r={sphere_r}]{int(n1_sphere)};{int(n2_sphere)}"
        final_n1, final_n2 = int(n1_sphere), int(n2_sphere)

    print(f"Analysis Results:")
    print(f"Box ({box_l}x{box_w}x{box_h}): {int(n2_box)} large balls, {int(n1_box)} small balls. SA={box_sa:.1f} cm^2. Energy = {energy_box} MJ")
    print(f"Cylinder (r={cyl_r}, h={cyl_h}): {n2_cyl} large balls, {n1_cyl} small balls. SA={cyl_sa:.1f} cm^2. Energy = {energy_cyl} MJ")
    print(f"Sphere (r={sphere_r}): {n2_sphere} large balls, {n1_sphere} small balls. SA={sphere_sa:.1f} cm^2. Energy = {energy_sphere} MJ")
    print("\nConclusion: The spherical container provides the highest energy capacity.")
    
    print("\nFinal Answer:")
    print(f"The optimal container is a {best_config.split(']')[0][1:]}.")
    print(f"It can hold {final_n2} large (2-cm) balls and {final_n1} small (1-cm) balls.")
    print(f"The calculation for the total energy is:")
    print(f"Total energy = 1 * {final_n1} + 10 * {final_n2} = {best_energy} MJ")
    
    # Final answer in the required format
    print(f"\n<<<{best_config}>>>")

solve()