import math
from collections import defaultdict

def solve_pizza_combinations():
    """
    Calculates the number of possible sets of 3, 4, and 5 pizza sizes
    that meet the specified franchising criteria.
    """
    
    # Use a defaultdict with 'set' to automatically handle unique diameters per slice area
    area_to_diameters = defaultdict(set)
    
    # Define constraints
    max_diameter = 20
    # A 20-inch pizza has an area of ~314 sq inches. With a min slice area of 5,
    # the theoretical max number of slices is ~62. We'll use 64 as a safe upper bound.
    max_slices = 64
    min_slice_area = 5.0
    max_slice_area = 20.0

    # Step 1 & 2: Iterate through all possibilities to find valid pizzas and group them by slice area
    for d in range(1, max_diameter + 1):
        for s in range(4, max_slices + 1, 4):
            # Area of a circle = pi * r^2 = pi * (d/2)^2 = (pi * d^2) / 4
            pizza_area = math.pi * (d ** 2) / 4
            slice_area = round(pizza_area / s, 1)
            
            # Check if the calculated slice area is within the allowed range
            if min_slice_area <= slice_area <= max_slice_area:
                area_to_diameters[slice_area].add(d)

    # Step 3: Count the combinations for sets of 3, 4, and 5
    total_combinations_3_sizes = 0
    total_combinations_4_sizes = 0
    total_combinations_5_sizes = 0

    for area in area_to_diameters:
        num_available_diameters = len(area_to_diameters[area])
        
        # Calculate combinations of 3 if enough diameters are available for this slice area
        if num_available_diameters >= 3:
            total_combinations_3_sizes += math.comb(num_available_diameters, 3)
            
        # Calculate combinations of 4
        if num_available_diameters >= 4:
            total_combinations_4_sizes += math.comb(num_available_diameters, 4)

        # Calculate combinations of 5
        if num_available_diameters >= 5:
            total_combinations_5_sizes += math.comb(num_available_diameters, 5)

    # Step 4: Sum the results to get the total number of possible sets
    total_possible_sets = total_combinations_3_sizes + total_combinations_4_sizes + total_combinations_5_sizes
    
    # Output the final counts as requested
    print(f"Number of possible 3-pizza sets: {total_combinations_3_sizes}")
    print(f"Number of possible 4-pizza sets: {total_combinations_4_sizes}")
    print(f"Number of possible 5-pizza sets: {total_combinations_5_sizes}")
    print(f"Total number of possible sets = {total_combinations_3_sizes} + {total_combinations_4_sizes} + {total_combinations_5_sizes} = {total_possible_sets}")

# Run the calculation and print the final answer
solve_pizza_combinations()