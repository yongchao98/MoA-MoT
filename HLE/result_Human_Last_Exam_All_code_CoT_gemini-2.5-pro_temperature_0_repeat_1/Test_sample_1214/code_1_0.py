import math

def solve_pioneer_storage():
    """
    Calculates the highest number of data points that can be collected and recorded
    by Pioneer's HDDs made from a 12x11cm rectangular material.
    """
    # Step 1: Define all the given parameters.
    
    # Material dimensions in cm
    rect_w = 12
    rect_h = 11

    # Disk 1 (D1) properties
    d1_diameter = 2
    d1_capacity_gb = 1

    # Disk 2 (D2) properties
    d2_diameter = 4
    d2_capacity_gb = 5

    # Step 2: Calculate the maximum possible storage capacity.
    # This function iterates through all grid-based packing combinations
    # to find the one that maximizes total storage capacity.
    max_capacity_gb = 0
    
    # We check both orientations of the material (12x11 and 11x12)
    for width, height in [(rect_w, rect_h), (rect_h, rect_w)]:
        max_d2_x = width // d2_diameter
        max_d2_y = height // d2_diameter
        
        # Iterate through all possible numbers of D2 disks in a rectangular block
        for i in range(max_d2_x + 1):
            for j in range(max_d2_y + 1):
                num_d2 = i * j
                
                # Calculate the area used by the block of D2 disks
                used_w = i * d2_diameter
                used_h = j * d2_diameter
                
                # The remaining area is an L-shape, which we can split into two rectangles
                # to fill with D1 disks.
                
                # Rectangle A: The strip next to the D2 block
                rem_w_A = width - used_w
                rem_h_A = height
                num_d1_A = (rem_w_A // d1_diameter) * (rem_h_A // d1_diameter)
                
                # Rectangle B: The strip below the D2 block
                rem_w_B = used_w
                rem_h_B = height - used_h
                num_d1_B = (rem_w_B // d1_diameter) * (rem_h_B // d1_diameter)
                
                num_d1 = num_d1_A + num_d1_B
                
                current_capacity = (num_d2 * d2_capacity_gb) + (num_d1 * d1_capacity_gb)
                
                if current_capacity > max_capacity_gb:
                    max_capacity_gb = current_capacity

    # Step 3: Estimate the size of a single data point.
    # The longest possible entry is "12 pm: extreme cold" (19 chars) plus a newline char (\n).
    # Assuming 1 character = 1 byte.
    data_point_size_bytes = 20

    # Step 4: Calculate the total number of data points.
    # Convert total capacity to bytes (1 GB = 1024^3 Bytes)
    gb_to_bytes = 1024**3
    max_capacity_bytes = max_capacity_gb * gb_to_bytes
    total_data_points = max_capacity_bytes // data_point_size_bytes

    # Step 5: Print the final equation and result.
    print(f"The maximum storage capacity achievable is {max_capacity_gb} GB.")
    print(f"The estimated size of one data point is {data_point_size_bytes} Bytes.")
    print("\nThe highest number of data points is calculated as follows:")
    print(f"Total Data Points = (Total Capacity in GB * Bytes per GB) / Bytes per Data Point")
    print(f"Total Data Points = ({max_capacity_gb} * {gb_to_bytes}) / {data_point_size_bytes} = {total_data_points}")

solve_pioneer_storage()
<<<1932735283>>>