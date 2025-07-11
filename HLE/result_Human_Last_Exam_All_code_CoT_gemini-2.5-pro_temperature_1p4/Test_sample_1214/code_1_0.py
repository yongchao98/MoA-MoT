import math

def solve_pioneer_storage():
    """
    Calculates the highest number of data points that can be collected and recorded by the Pioneer probe.
    """
    # Step 1: Define initial parameters
    material_w_cm = 12
    material_h_cm = 11

    # D1 disk properties
    d1_radius_cm = 1
    d1_diameter_cm = d1_radius_cm * 2
    d1_capacity_gb = 1

    # D2 disk properties
    d2_radius_cm = 2
    d2_diameter_cm = d2_radius_cm * 2
    d2_capacity_gb = 5

    # Step 2: Determine the optimal disk layout to maximize capacity.
    # We prioritize the larger, more capacious D2 disks.
    # We place a grid of D2 disks and then fill the remaining area with D1 disks.

    # Fit D2 disks (4cm diameter) onto the 12x11 cm material
    d2_fit_along_w = material_w_cm // d2_diameter_cm
    d2_fit_along_h = material_h_cm // d2_diameter_cm
    num_d2_disks = d2_fit_along_w * d2_fit_along_h

    # Calculate remaining area for D1 disks
    # The D2 disks occupy a 12cm x 8cm area (3x2 grid)
    area_used_by_d2_h = d2_fit_along_h * d2_diameter_cm
    remaining_material_w = material_w_cm
    remaining_material_h = material_h_cm - area_used_by_d2_h

    # Fit D1 disks (2cm diameter) into the remaining 12x3 cm area
    if remaining_material_h > 0:
      d1_fit_along_rem_w = remaining_material_w // d1_diameter_cm
      d1_fit_along_rem_h = remaining_material_h // d1_diameter_cm
      num_d1_disks = d1_fit_along_rem_w * d1_fit_along_rem_h
    else:
      num_d1_disks = 0


    # Step 3: Calculate total storage capacity
    total_capacity_gb = (num_d1_disks * d1_capacity_gb) + (num_d2_disks * d2_capacity_gb)
    
    # Step 4: Estimate data point size
    # Example data point: "16 pm: extreme hot" + newline
    # Time (max 5 chars) + ": " (2 chars) + Category (max 12 for "extreme cold") + newline (1 char)
    # Assuming 1 character = 1 byte (ASCII/UTF-8 for simple characters)
    size_per_datapoint_bytes = 5 + 2 + 12 + 1
    
    # Step 5: Convert total capacity to bytes and calculate total data points
    bytes_per_gb = 1024**3
    total_capacity_bytes = total_capacity_gb * bytes_per_gb
    
    total_data_points = total_capacity_bytes // size_per_datapoint_bytes

    # --- Print the results and the final equation ---
    print("--- Solving Pioneer's Storage Problem ---")
    print(f"\n1. Optimal Disk Configuration from a {material_w_cm}x{material_h_cm} cm sheet:")
    print(f"   - Number of D2 disks (5GB): {num_d2_disks}")
    print(f"   - Number of D1 disks (1GB): {num_d1_disks}")
    
    print("\n2. Total Storage Capacity Calculation:")
    print(f"   - Capacity from D1 disks: {num_d1_disks} * {d1_capacity_gb} GB = {num_d1_disks * d1_capacity_gb} GB")
    print(f"   - Capacity from D2 disks: {num_d2_disks} * {d2_capacity_gb} GB = {num_d2_disks * d2_capacity_gb} GB")
    print(f"   - Total Capacity: {total_capacity_gb} GB")

    print("\n3. Final Calculation for Data Points:")
    print(f"   - Size of one data point: {size_per_datapoint_bytes} bytes")
    print(f"   - Bytes per Gigabyte: {bytes_per_gb}")
    
    print("\nThe final equation to calculate the total number of data points is:")
    print(f"({num_d2_disks} * {d2_capacity_gb} + {num_d1_disks} * {d1_capacity_gb}) * {bytes_per_gb} // {size_per_datapoint_bytes} = {total_data_points}")

if __name__ == "__main__":
    solve_pioneer_storage()
    # Final answer in the required format
    print("\n<<<1932735283>>>")