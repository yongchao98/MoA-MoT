import math

def solve_pioneer_storage():
    """
    Calculates the maximum number of data points (interpreted as GB of storage)
    that can be collected by Pioneer.
    """
    # Material dimensions in cm
    mat_w = 12
    mat_h = 11

    # Disk 1 (D1) properties
    d1_diameter = 2
    d1_capacity = 1  # in GB

    # Disk 2 (D2) properties
    d2_diameter = 4
    d2_capacity = 5  # in GB

    # --- Calculation ---
    # The strategy is to prioritize the higher-capacity D2 disks.
    # We will fit as many D2 disks as possible and use the remainder for D1 disks.

    # Step 1: Calculate the maximum number of D2 disks (4cm diameter) from the 12x11 cm material.
    num_d2_along_w = mat_w // d2_diameter
    num_d2_along_h = mat_h // d2_diameter
    num_d2 = num_d2_along_w * num_d2_along_h

    # Step 2: Calculate the dimensions of the material remaining after cutting the D2 disks.
    # The D2 disks occupy a (3*4) x (2*4) = 12cm x 8cm area.
    used_h_by_d2 = num_d2_along_h * d2_diameter
    rem_mat_w = mat_w
    rem_mat_h = mat_h - used_h_by_d2

    # Step 3: Calculate the maximum number of D1 disks (2cm diameter) from the remaining material.
    # The remaining piece is 12cm x 3cm.
    num_d1_along_rem_w = rem_mat_w // d1_diameter
    num_d1_along_rem_h = rem_mat_h // d1_diameter
    num_d1 = num_d1_along_rem_w * num_d1_along_rem_h

    # Step 4: Calculate the total storage capacity.
    total_capacity = (num_d2 * d2_capacity) + (num_d1 * d1_capacity)
    
    # As requested, printing the final equation with each number.
    # The "number of data points" is interpreted as the total storage in GB.
    print("Based on the optimal cutting strategy:")
    print(f"Number of D2 disks (5GB each): {num_d2}")
    print(f"Number of D1 disks (1GB each): {num_d1}")
    print("\nThe final equation for the total number of data points (in GB) is:")
    print(f"({num_d2} * {d2_capacity}) + ({num_d1} * {d1_capacity}) = {total_capacity}")

solve_pioneer_storage()