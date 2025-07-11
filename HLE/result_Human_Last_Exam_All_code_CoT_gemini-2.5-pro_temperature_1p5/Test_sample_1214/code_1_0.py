import math

def calculate_max_data_points():
    """
    Calculates the highest number of data points that can be collected and recorded based on
    the storage capacity that can be manufactured from a given material.
    """
    # --- Step 1: Define initial parameters ---
    rect_length_cm = 12
    rect_width_cm = 11

    # Disk type 1 specifications
    d1_radius_cm = 1
    d1_capacity_gb = 1

    # Disk type 2 specifications
    d2_radius_cm = 2
    d2_capacity_gb = 5

    # Data point size estimation (e.g., "12 am: extreme cold\n" is 21 chars -> 21 bytes)
    data_point_size_bytes = 21
    
    # Conversion factor for Gigabytes to Bytes (1 GiB = 2^30 bytes)
    BYTES_PER_GB = 1024 * 1024 * 1024

    print("Step 1: Calculating the number of disks that can be made from the 12x11 cm material.")
    print("-" * 70)

    # --- Step 2: Calculate max disks of each type using grid packing ---
    d1_diameter_cm = 2 * d1_radius_cm
    d2_diameter_cm = 2 * d2_radius_cm

    # For D1 disks (diameter 2cm)
    num_d1_orient1 = math.floor(rect_length_cm / d1_diameter_cm) * math.floor(rect_width_cm / d1_diameter_cm)
    num_d1_orient2 = math.floor(rect_width_cm / d1_diameter_cm) * math.floor(rect_length_cm / d1_diameter_cm)
    max_d1_disks = max(num_d1_orient1, num_d1_orient2)
    
    print(f"For D1 disks (radius {d1_radius_cm}cm):")
    print(f"  - We can fit {math.floor(rect_length_cm/d1_diameter_cm)} disks along the 12cm side and {math.floor(rect_width_cm/d1_diameter_cm)} along the 11cm side.")
    print(f"  - Total D1 disks = {max_d1_disks}")

    # For D2 disks (diameter 4cm)
    num_d2_orient1 = math.floor(rect_length_cm / d2_diameter_cm) * math.floor(rect_width_cm / d2_diameter_cm)
    num_d2_orient2 = math.floor(rect_width_cm / d2_diameter_cm) * math.floor(rect_length_cm / d2_diameter_cm)
    max_d2_disks = max(num_d2_orient1, num_d2_orient2)

    print(f"\nFor D2 disks (radius {d2_radius_cm}cm):")
    print(f"  - We can fit {math.floor(rect_length_cm/d2_diameter_cm)} disks along the 12cm side and {math.floor(rect_width_cm/d2_diameter_cm)} along the 11cm side.")
    print(f"  - Total D2 disks = {max_d2_disks}")
    
    print("\nStep 2: Calculating total storage capacity for each scenario.")
    print("-" * 70)
    
    # --- Step 3: Calculate total capacity for each scenario ---
    total_capacity_d1_gb = max_d1_disks * d1_capacity_gb
    total_capacity_d2_gb = max_d2_disks * d2_capacity_gb

    print(f"Total capacity with D1 disks: {max_d1_disks} disks * {d1_capacity_gb} GB/disk = {total_capacity_d1_gb} GB")
    print(f"Total capacity with D2 disks: {max_d2_disks} disks * {d2_capacity_gb} GB/disk = {total_capacity_d2_gb} GB")
    
    # --- Step 4: Determine the maximum capacity ---
    max_capacity_gb = max(total_capacity_d1_gb, total_capacity_d2_gb)
    print(f"\nThe maximum possible storage capacity is {max_capacity_gb} GB.")

    print("\nStep 3: Calculating the total number of data points.")
    print("-" * 70)

    # --- Step 5: Calculate total number of data points ---
    max_capacity_bytes = max_capacity_gb * BYTES_PER_GB
    total_data_points = math.floor(max_capacity_bytes / data_point_size_bytes)
    
    print(f"Assumptions:")
    print(f"  - Size of one data point (e.g., '12 am: extreme cold\\n'): {data_point_size_bytes} bytes")
    print(f"  - 1 GB = {BYTES_PER_GB} bytes\n")
    
    print("Final Calculation:")
    print(f"Highest number of data points = (Max Capacity in Bytes) / (Size per Data Point)")
    final_equation_str = (
        f"                               = ({max_capacity_gb} GB * {BYTES_PER_GB} bytes/GB) / {data_point_size_bytes} bytes"
    )
    print(final_equation_str)
    print(f"                               = {max_capacity_bytes} / {data_point_size_bytes}")
    print(f"                               = {total_data_points}")
    
    return total_data_points

if __name__ == "__main__":
    final_answer = calculate_max_data_points()
    print(f"\n<<<The highest number of data points that can be collected and recorded is {final_answer}.>>>")
    # This format is just for the final answer extraction.
    print(f"\n<<<{final_answer}>>>")