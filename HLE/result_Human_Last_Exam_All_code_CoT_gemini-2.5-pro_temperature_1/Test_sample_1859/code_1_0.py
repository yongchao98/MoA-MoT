import math

def find_pizza_combinations():
    """
    Calculates the number of possible sets of three, four, and five pizza sizes
    based on a series of constraints on diameter, slices, and slice area.
    """
    # Define the problem constraints
    max_diameter = 20
    min_slice_area = 5.0
    max_slice_area = 20.0

    # Step 1 & 2: Find all valid individual pizza configurations and group them by slice area.
    # This dictionary will store {slice_area: [list_of_diameters]}
    valid_pizzas_by_area = {}

    # Iterate through all possible integer diameters from 1 to 20 inches
    for d in range(1, max_diameter + 1):
        total_area = math.pi * (d / 2)**2

        # To avoid division by zero later and handle impossible pizzas
        if total_area < min_slice_area:
            continue
            
        # Determine the valid range for the number of slices (s) based on area constraints
        # min_slice_area <= total_area / s  => s <= total_area / min_slice_area
        # max_slice_area >= total_area / s  => s >= total_area / max_slice_area
        min_s = math.ceil(total_area / max_slice_area)
        max_s = math.floor(total_area / min_slice_area)
        
        # Find the first valid slice count (must be a multiple of 4)
        start_s = ((min_s + 3) // 4) * 4
        
        # Iterate through possible slice counts in steps of 4
        for s in range(start_s, max_s + 1, 4):
            if s == 0:
                continue
                
            # Calculate slice area, rounded to one decimal place, as required
            slice_area = round(total_area / s, 1)
            
            # Check if the calculated slice area is within the allowed range
            if min_slice_area <= slice_area <= max_slice_area:
                # Group the diameter by its resulting slice area
                if slice_area not in valid_pizzas_by_area:
                    valid_pizzas_by_area[slice_area] = []
                
                # Add the diameter to the list, ensuring it's unique for that area
                if d not in valid_pizzas_by_area[slice_area]:
                    valid_pizzas_by_area[slice_area].append(d)

    # Step 3 & 4: Calculate the number of combinations for sets of 3, 4, and 5
    combinations_of_3 = 0
    combinations_of_4 = 0
    combinations_of_5 = 0

    for area in valid_pizzas_by_area:
        diameters = valid_pizzas_by_area[area]
        n = len(diameters)
        
        # Calculate combinations for sets of three
        if n >= 3:
            combinations_of_3 += math.comb(n, 3)
            
        # Calculate combinations for sets of four
        if n >= 4:
            combinations_of_4 += math.comb(n, 4)

        # Calculate combinations for sets of five
        if n >= 5:
            combinations_of_5 += math.comb(n, 5)

    # Step 5: Sum and print the results
    total_combinations = combinations_of_3 + combinations_of_4 + combinations_of_5
    
    print(f"Possible sets of three pizza sizes: {combinations_of_3}")
    print(f"Possible sets of four pizza sizes: {combinations_of_4}")
    print(f"Possible sets of five pizza sizes: {combinations_of_5}")
    print(f"Total possible sets = {combinations_of_3} + {combinations_of_4} + {combinations_of_5} = {total_combinations}")
    
    return total_combinations

# Run the function to get the answer
final_answer = find_pizza_combinations()
# The final answer is wrapped for automated reading.
# print(f"<<<{final_answer}>>>")