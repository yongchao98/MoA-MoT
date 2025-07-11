import math

def solve():
    """
    Calculates the maximum number of data points that can be stored on HDDs
    manufactured from a given sheet of material.
    """
    # Define material and disk specifications
    sheet_width_cm = 12
    sheet_height_cm = 11

    d1_radius_cm = 1
    d1_capacity_gb = 1
    d1_diameter_cm = d1_radius_cm * 2

    d2_radius_cm = 2
    d2_capacity_gb = 5
    d2_diameter_cm = d2_radius_cm * 2

    # --- Step 1: Estimate the size of a single data point ---
    # A data point like "16 pm: extreme hot\n" is stored.
    # Let's estimate its size in bytes:
    # - Time ("16 pm"): ~8 bytes for a flexible format
    # - Separator (": "): 2 bytes
    # - Category (longest is "extreme cold"): ~16 bytes for string and null terminator
    # - Newline character: 1 byte
    # This gives a total of 8 + 2 + 16 + 1 = 27 bytes.
    # We will round up to 32 bytes to account for file system alignment and overhead.
    datapoint_size_bytes = 32

    # --- Step 2: Determine the optimal disk arrangement to maximize storage ---
    # D2 disks are more capacity-dense (5GB for a 4cm diameter circle) than D1 disks (1GB for a 2cm one).
    # So, we should maximize the number of D2 disks first.
    # We use a simple grid packing model.

    # Calculate how many D2 disks can fit on the 12x11cm sheet.
    num_d2_cols = sheet_width_cm // d2_diameter_cm
    num_d2_rows = sheet_height_cm // d2_diameter_cm
    num_d2_disks = num_d2_cols * num_d2_rows

    # The D2 disks will occupy an area of (num_d2_cols * d2_diameter_cm) by (num_d2_rows * d2_diameter_cm).
    # This is (3 * 4cm) x (2 * 4cm) = 12cm x 8cm.
    # The remaining material is a rectangular strip.
    remaining_width_cm = sheet_width_cm
    remaining_height_cm = sheet_height_cm - (num_d2_rows * d2_diameter_cm)

    # Now, calculate how many D1 disks can be cut from the remaining strip (12cm x 3cm).
    if remaining_width_cm >= d1_diameter_cm and remaining_height_cm >= d1_diameter_cm:
        num_d1_disks = (remaining_width_cm // d1_diameter_cm) * (remaining_height_cm // d1_diameter_cm)
    else:
        num_d1_disks = 0
    
    # --- Step 3: Calculate total storage capacity ---
    total_capacity_gb = (num_d2_disks * d2_capacity_gb) + (num_d1_disks * d1_capacity_gb)

    # Conversion factor from Gigabytes to Bytes (1 GB = 1024^3 Bytes)
    bytes_per_gb = 1024**3
    total_capacity_bytes = total_capacity_gb * bytes_per_gb

    # --- Step 4: Calculate total number of data points ---
    max_data_points = total_capacity_bytes // datapoint_size_bytes

    # --- Print the final equation and result ---
    print("This problem is solved by maximizing storage and then dividing by the size of one data point.")
    print(f"1. Maximize Storage: Cut {num_d2_disks} D2 disks and {num_d1_disks} D1 disks from the sheet.")
    print(f"2. Total Capacity (GB) = ({num_d2_disks} * {d2_capacity_gb} GB) + ({num_d1_disks} * {d1_capacity_gb} GB) = {total_capacity_gb} GB")
    print(f"3. Estimate Data Point Size: {datapoint_size_bytes} bytes per data point.")
    print("\nFinal calculation combines these steps:")
    # Print the equation with all the numbers, as requested
    print(f"(({num_d2_disks} * {d2_capacity_gb}) + ({num_d1_disks} * {d1_capacity_gb})) * {bytes_per_gb} / {datapoint_size_bytes} = {max_data_points}")

solve()
<<<1207959552>>>