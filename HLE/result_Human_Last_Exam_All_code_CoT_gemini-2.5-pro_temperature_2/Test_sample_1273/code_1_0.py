import math

def solve_billet_problem():
    """
    Calculates the maximum value that can be obtained by cutting a steel billet
    into balls and cubes of specified dimensions and prices.
    """

    # --- Problem Parameters ---
    
    # Billet dimensions in cm
    billet_length = 8
    billet_width = 8
    billet_height = 4

    # Product B2: Ball of 2cm radius
    b2_radius = 2
    b2_price = 150

    # Product T1: Cube of 0.8cm side
    t1_side = 0.8
    t1_price = 1

    # Product B1: Ball of 1cm diameter (used for comparison)
    b1_diameter = 1
    b1_radius = b1_diameter / 2
    b1_price = 1

    # --- Step 1: Maximize the number of B2 balls ---
    
    # Each B2 ball (2cm radius) has a 4cm diameter, requiring a 4x4x4 cm bounding box.
    b2_diameter = 2 * b2_radius
    num_b2 = math.floor(billet_length / b2_diameter) * \
             math.floor(billet_width / b2_diameter) * \
             math.floor(billet_height / b2_diameter)

    value_from_b2 = num_b2 * b2_price

    # --- Step 2: Calculate the volume of leftover material (scrap) ---
    
    billet_volume = billet_length * billet_width * billet_height
    # Volume of a single B2 sphere: (4/3) * pi * r^3
    b2_sphere_volume = (4/3) * math.pi * (b2_radius ** 3)
    total_volume_of_b2_spheres = num_b2 * b2_sphere_volume
    scrap_volume = billet_volume - total_volume_of_b2_spheres

    # --- Step 3: Maximize value from the scrap volume ---
    
    # Compare volumes of T1 and B1 to see which is a better use of scrap.
    # Price is the same, so the smaller volume is better.
    t1_cube_volume = t1_side ** 3  # 0.8^3 = 0.512
    b1_sphere_volume = (4/3) * math.pi * (b1_radius ** 3) # ~0.5236
    # T1 cubes are smaller, so we choose them.
    
    num_t1_from_scrap = math.floor(scrap_volume / t1_cube_volume)
    value_from_t1 = num_t1_from_scrap * t1_price

    # --- Step 4: Calculate the total maximum value ---

    total_value = value_from_b2 + value_from_t1

    # --- Step 5: Output the final equation as requested ---
    
    print("Based on the calculations, the final equation for the maximum value is:")
    print(f"Total Value = ({num_b2} * {b2_price}) + ({num_t1_from_scrap} * {t1_price}) = {total_value}")


solve_billet_problem()