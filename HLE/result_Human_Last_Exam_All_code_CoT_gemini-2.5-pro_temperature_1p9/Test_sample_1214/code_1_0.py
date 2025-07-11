import math

def solve_pioneer_storage():
    """
    Calculates the maximum number of data points that can be stored on HDDs
    manufactured from a given sheet of material.
    """

    # Step 1: Define constants and estimate data point size
    
    # Material dimensions in cm
    sheet_w = 12
    sheet_h = 11

    # Disk properties
    # D1: 1cm radius -> 2cm diameter, requires a 2x2 cm square
    d1_diameter = 2
    d1_capacity_gb = 1

    # D2: 2cm radius -> 4cm diameter, requires a 4x4 cm square
    d2_diameter = 4
    d2_capacity_gb = 5

    # Storage and data point size definitions
    # Longest data point string: "12 am: extreme cold" (19 chars) + 1 for newline = 20 bytes
    data_point_size_bytes = 20
    # Using manufacturer's definition for Gigabyte (1 billion bytes)
    bytes_per_gb = 10**9

    # Step 2: Solve the packing problem to maximize storage capacity
    max_capacity_gb = 0
    
    # Strategy A: Only D1 disks (2x2 cm)
    num_d1_only = math.floor(sheet_w / d1_diameter) * math.floor(sheet_h / d1_diameter)
    capacity_a = num_d1_only * d1_capacity_gb
    
    # Strategy B: Only D2 disks (4x4 cm)
    num_d2_only = math.floor(sheet_w / d2_diameter) * math.floor(sheet_h / d2_diameter)
    capacity_b = num_d2_only * d2_capacity_gb

    # Strategy C: Mixed - Place D2 first, then fill remainder with D1
    # Place a 3x2 grid of D2 disks (4x4 cm each)
    num_d2_mixed = math.floor(sheet_w / d2_diameter) * math.floor(sheet_h / d2_diameter) # 3 * 2 = 6 D2 disks
    # This configuration uses a 12cm x 8cm area
    area_used_w = num_d2_mixed // math.floor(sheet_h/d2_diameter) * d2_diameter # 3 * 4 = 12
    area_used_h = math.floor(sheet_h / d2_diameter) * d2_diameter # 2 * 4 = 8
    
    # The remaining area is a 12cm x 3cm strip
    rem_w = sheet_w 
    rem_h = sheet_h - area_used_h
    
    # Fit D1 disks into the remaining strip
    num_d1_mixed = math.floor(rem_w / d1_diameter) * math.floor(rem_h / d1_diameter)
    
    capacity_c = (num_d2_mixed * d2_capacity_gb) + (num_d1_mixed * d1_capacity_gb)
    
    # The optimal capacity is the maximum of the calculated strategies
    max_capacity_gb = max(capacity_a, capacity_b, capacity_c)

    # Step 3 & 4: Calculate total bytes and the final number of data points
    total_capacity_bytes = max_capacity_gb * bytes_per_gb
    total_data_points = total_capacity_bytes // data_point_size_bytes

    # Print the final result including the equation
    print("This plan finds the optimal mix of disks to maximize storage and calculates the total data points.")
    print(f"Optimal Disk Configuration: {num_d2_mixed} D2 disks and {num_d1_mixed} D1 disks.")
    print(f"Maximum Storage Capacity: {max_capacity_gb} GB.")
    print("\nFinal Calculation:")
    print(f"Highest number of data points = ({max_capacity_gb} GB * {bytes_per_gb:,} bytes/GB) / {data_point_size_bytes} bytes/datapoint")
    print(f"Highest number of data points = {total_data_points:,}")
    
    # Final answer in the required format
    print(f"\n<<<{total_data_points}>>>")

solve_pioneer_storage()