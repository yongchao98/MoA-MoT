import math

def solve_pioneer_storage():
    """
    Calculates the maximum number of data points Pioneer can record based on
    HDD capacity derived from a given material sheet.
    """
    # --- Part 1: Define constants and problem parameters ---
    mat_w = 12  # cm
    mat_h = 11  # cm

    d1_radius = 1  # cm
    d1_capacity_gb = 1  # GB

    d2_radius = 2  # cm
    d2_capacity_gb = 5  # GB

    # Size of one data point is based on the longest example: "12 am: extreme cold"
    # Plus a colon, a space, and a newline character for separation.
    # "12 am" (5) + ": " (2) + "extreme cold" (12) + "\n" (1) = 20 bytes
    data_point_size_bytes = 20

    # Standard conversion for Gigabytes to Bytes
    bytes_per_gb = 1024**3

    # Calculate the side length of the square material needed for each disk (diameter)
    d1_side = d1_radius * 2
    d2_side = d2_radius * 2

    # --- Part 2: Calculate the best disk layout for maximum capacity ---
    # We test two orientations of the material sheet (12x11 and 11x12)
    # to find the optimal packing, prioritizing the larger D2 disks.

    # Layout 1: Material is 12cm x 11cm
    # Fit D2 disks (4x4 cm squares)
    num_d2_w1 = mat_w // d2_side  # 12 // 4 = 3
    num_d2_h1 = mat_h // d2_side  # 11 // 4 = 2
    num_d2_l1 = num_d2_w1 * num_d2_h1

    # The remaining area is a 12cm x 3cm strip
    rem_area_w1 = mat_w
    rem_area_h1 = mat_h - (num_d2_h1 * d2_side)
    # Fit D1 disks (2x2 cm squares) into the remaining strip
    num_d1_l1 = (rem_area_w1 // d1_side) * (rem_area_h1 // d1_side)
    cap_l1 = (num_d2_l1 * d2_capacity_gb) + (num_d1_l1 * d1_capacity_gb)

    # Layout 2: Material is 11cm x 12cm
    # Fit D2 disks (4x4 cm squares)
    num_d2_w2 = mat_h // d2_side # 11 // 4 = 2
    num_d2_h2 = mat_w // d2_side # 12 // 4 = 3
    num_d2_l2 = num_d2_w2 * num_d2_h2

    # The remaining area is a 3cm x 12cm strip
    rem_area_w2 = mat_h - (num_d2_w2 * d2_side)
    rem_area_h2 = mat_w
    # Fit D1 disks (2x2 cm squares) into the remaining strip
    num_d1_l2 = (rem_area_w2 // d1_side) * (rem_area_h2 // d1_side)
    cap_l2 = (num_d2_l2 * d2_capacity_gb) + (num_d1_l2 * d1_capacity_gb)

    # Determine the best configuration
    if cap_l1 > cap_l2:
        max_capacity_gb = cap_l1
        num_d2_best = num_d2_l1
        num_d1_best = num_d1_l1
    else: # cap_l2 >= cap_l1
        max_capacity_gb = cap_l2
        num_d2_best = num_d2_l2
        num_d1_best = num_d1_l1

    # --- Part 3: Calculate the total number of data points ---
    total_capacity_bytes = max_capacity_gb * bytes_per_gb
    total_data_points = total_capacity_bytes // data_point_size_bytes

    # --- Part 4: Print the detailed results ---
    print("### Step-by-step Calculation ###\n")
    print("--- 1. Maximizing Storage Capacity ---")
    print(f"The manufacturing material is a {mat_w}cm x {mat_h}cm rectangle.")
    print(f"To maximize capacity, we fit the largest disks (D2, {d2_capacity_gb}GB from {d2_side}x{d2_side}cm squares) first.")
    print(f"The best configuration found uses {num_d2_best} D2 disks and {num_d1_best} D1 disks.")

    d2_cap_total = num_d2_best * d2_capacity_gb
    d1_cap_total = num_d1_best * d1_capacity_gb
    print(f"Capacity from D2 disks: {num_d2_best} disks * {d2_capacity_gb} GB/disk = {d2_cap_total} GB")
    print(f"Capacity from D1 disks: {num_d1_best} disks * {d1_capacity_gb} GB/disk = {d1_cap_total} GB")
    print(f"Maximum Total Capacity = {d2_cap_total} + {d1_cap_total} = {max_capacity_gb} GB\n")

    print("--- 2. Calculating Total Data Points ---")
    print(f"Total capacity in Bytes = {max_capacity_gb} GB * {bytes_per_gb:,} Bytes/GB = {total_capacity_bytes:,} Bytes")
    print(f"Size of one data point = {data_point_size_bytes} Bytes\n")

    print("--- 3. Final Equation ---")
    print(f"Total data points = Total Capacity (Bytes) / Size per data point (Bytes)")
    print(f"Equation: {total_capacity_bytes} / {data_point_size_bytes}")
    print(f"Result = {total_data_points:,} data points")

    print(f"\n<<<{total_data_points}>>>")

solve_pioneer_storage()