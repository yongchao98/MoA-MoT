import math

def calculate_max_data_points():
    """
    Calculates the highest number of data points that can be collected and recorded.
    """
    # Step 1: Define material and disk properties
    rect_w = 12  # cm
    rect_h = 11  # cm

    d1_radius = 1  # cm
    d1_capacity_gb = 1  # GB

    d2_radius = 2  # cm
    d2_capacity_gb = 5  # GB

    d1_diameter = d1_radius * 2
    d2_diameter = d2_radius * 2

    # Step 2: Calculate the optimal number of disks to maximize capacity.
    # Strategy: Pack the larger, more capacity-dense D2 disks first, then fill the rest with D1.
    
    # Pack D2 disks (4cm diameter) in a 12x11 cm area
    num_d2_along_12 = math.floor(rect_w / d2_diameter)
    num_d2_along_11 = math.floor(rect_h / d2_diameter)
    num_d2_disks = num_d2_along_12 * num_d2_along_11

    # Calculate the remaining area after packing D2 disks
    area_used_by_d2_w = num_d2_along_12 * d2_diameter
    area_used_by_d2_h = num_d2_along_11 * d2_diameter
    
    # The remaining area is a rectangle of size (rect_w) x (rect_h - area_used_by_d2_h)
    rem_area_w = rect_w
    rem_area_h = rect_h - area_used_by_d2_h

    # Pack D1 disks (2cm diameter) in the remaining area
    num_d1_along_rem_w = math.floor(rem_area_w / d1_diameter)
    num_d1_along_rem_h = math.floor(rem_area_h / d1_diameter)
    num_d1_disks = num_d1_along_rem_w * num_d1_along_rem_h

    # Step 3: Calculate the total storage capacity
    total_capacity_gb = (num_d2_disks * d2_capacity_gb) + (num_d1_disks * d1_capacity_gb)

    # Step 4: Calculate the number of data points
    # Assume a single data point record size is 32 bytes.
    data_point_size_bytes = 32
    bytes_per_gb = 1024**3
    
    total_capacity_bytes = total_capacity_gb * bytes_per_gb
    
    total_data_points = math.floor(total_capacity_bytes / data_point_size_bytes)

    # Print the results step-by-step
    print("--- Disk Manufacturing Calculation ---")
    print(f"From a {rect_w}x{rect_h} cm sheet:")
    print(f"Number of D2 disks (4cm diameter): {num_d2_along_12} x {num_d2_along_11} = {num_d2_disks}")
    print(f"Number of D1 disks (2cm diameter) in remaining space: {num_d1_along_rem_w} x {num_d1_along_rem_h} = {num_d1_disks}")
    print("\n--- Capacity Calculation ---")
    print(f"Total Capacity = ({num_d2_disks} D2 disks * {d2_capacity_gb} GB) + ({num_d1_disks} D1 disks * {d1_capacity_gb} GB)")
    print(f"Total Capacity = {num_d2_disks * d2_capacity_gb} GB + {num_d1_disks * d1_capacity_gb} GB = {total_capacity_gb} GB")
    
    print("\n--- Data Point Calculation ---")
    print(f"Assuming 1 data point = {data_point_size_bytes} bytes")
    print(f"Total Capacity in Bytes = {total_capacity_gb} GB * {bytes_per_gb} B/GB = {total_capacity_bytes} Bytes")
    print("\n--- Final Equation ---")
    print(f"Highest number of data points = Total Capacity (Bytes) / Size per data point (Bytes)")
    print(f"Highest number of data points = {total_capacity_bytes} / {data_point_size_bytes} = {total_data_points}")

calculate_max_data_points()
<<<1207959552>>>