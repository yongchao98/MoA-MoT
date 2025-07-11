import math

def solve_pizza_combinations():
    """
    Calculates the number of possible sets of 3, 4, and 5 pizza sizes
    based on a given set of constraints.
    """

    # A dictionary to store {rounded_slice_area: {set_of_diameters}}
    area_to_diameters = {}
    
    # Define the constraints from the problem description
    min_diameter = 1
    max_diameter = 20
    min_area = 5.0
    max_area = 20.0
    
    # Step 1 & 2: Find all valid pizza configurations and group them by slice area.
    # Iterate through all possible integer diameters.
    for diameter in range(min_diameter, max_diameter + 1):
        # Iterate through possible numbers of slices (multiples of 4).
        # An upper limit of 64 slices is sufficient, as larger values
        # with the max diameter would result in a slice area below the minimum.
        for num_slices in range(4, 65, 4):
            # Calculate the area of a single slice
            radius = diameter / 2.0
            total_area = math.pi * (radius ** 2)
            slice_area = total_area / num_slices
            
            # Round the area to one decimal place as required
            rounded_slice_area = round(slice_area, 1)
            
            # Check if the rounded area is within the valid range
            if min_area <= rounded_slice_area <= max_area:
                # If this area is new, initialize a set for it
                if rounded_slice_area not in area_to_diameters:
                    area_to_diameters[rounded_slice_area] = set()
                # Add the current diameter to the set for this slice area
                area_to_diameters[rounded_slice_area].add(diameter)

    # Step 3: Calculate the number of combinations for each group.
    total_comb_3 = 0
    total_comb_4 = 0
    total_comb_5 = 0
    
    # Iterate through the sets of diameters for each valid slice area
    for diameters in area_to_diameters.values():
        num_valid_diameters = len(diameters)
        
        # Calculate combinations of 3 pizza sizes
        if num_valid_diameters >= 3:
            total_comb_3 += math.comb(num_valid_diameters, 3)
            
        # Calculate combinations of 4 pizza sizes
        if num_valid_diameters >= 4:
            total_comb_4 += math.comb(num_valid_diameters, 4)
            
        # Calculate combinations of 5 pizza sizes
        if num_valid_diameters >= 5:
            total_comb_5 += math.comb(num_valid_diameters, 5)

    # Step 4: Sum the results and print the final output.
    total_combinations = total_comb_3 + total_comb_4 + total_comb_5
    
    print(f"Possible sets of 3 pizza sizes: {total_comb_3}")
    print(f"Possible sets of 4 pizza sizes: {total_comb_4}")
    print(f"Possible sets of 5 pizza sizes: {total_comb_5}")
    print("-" * 35)
    print(f"Total possible sets = {total_comb_3} + {total_comb_4} + {total_comb_5} = {total_combinations}")

# Run the function to get the answer
solve_pizza_combinations()
<<<870>>>