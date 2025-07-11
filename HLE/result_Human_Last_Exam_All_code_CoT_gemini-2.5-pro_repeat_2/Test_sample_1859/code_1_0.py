import math

def solve_pizza_combinations():
    """
    Calculates the number of possible sets of 3, 4, and 5 pizza sizes
    that meet the specified franchising criteria.
    """
    # Step 1 & 2: Identify possible pizzas and group them by slice area.
    # The dictionary will map a rounded_slice_area to a set of valid diameters.
    area_to_diameters = {}
    
    # Maximum possible number of slices. If d=20 (max) and slice_area=5 (min),
    # 5 = (pi * 20^2) / (4 * s) => s = (pi * 400) / 20 = 20 * pi ~= 62.8.
    # So, checking slices up to 60 (the nearest multiple of 4) is sufficient.
    max_slices = 60
    
    # Iterate through all allowed diameters (1-20 inches).
    for diameter in range(1, 21):
        # Iterate through all allowed slice numbers (multiples of 4).
        for slices in range(4, max_slices + 1, 4):
            # Calculate the actual area of a single slice.
            # Area = (pi * r^2) / slices = (pi * (d/2)^2) / slices
            actual_slice_area = (math.pi * (diameter ** 2)) / (4 * slices)
            
            # Round the area to one decimal place, as per the rules.
            rounded_slice_area = round(actual_slice_area, 1)
            
            # Check if the rounded slice area is within the allowed range [5, 20].
            if 5.0 <= rounded_slice_area <= 20.0:
                # If this area hasn't been seen before, initialize a set for it.
                if rounded_slice_area not in area_to_diameters:
                    area_to_diameters[rounded_slice_area] = set()
                # Add the current diameter to the set for this slice area.
                area_to_diameters[rounded_slice_area].add(diameter)

    # Step 3: Calculate combinations for each group of diameters.
    
    # Initialize counters for each set size.
    count_3_pizza_sets = 0
    count_4_pizza_sets = 0
    count_5_pizza_sets = 0
    
    # Helper function for combinations (n choose k).
    def combinations(n, k):
        if k < 0 or k > n:
            return 0
        return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

    # Iterate through the sets of diameters found for each common slice area.
    for diameters in area_to_diameters.values():
        num_diameters = len(diameters)
        
        # Calculate and add the number of possible sets of 3, 4, and 5.
        count_3_pizza_sets += combinations(num_diameters, 3)
        count_4_pizza_sets += combinations(num_diameters, 4)
        count_5_pizza_sets += combinations(num_diameters, 5)
        
    # Step 4: Sum the results and print the output.
    total_sets = count_3_pizza_sets + count_4_pizza_sets + count_5_pizza_sets
    
    print(f"Number of possible 3-pizza sets: {count_3_pizza_sets}")
    print(f"Number of possible 4-pizza sets: {count_4_pizza_sets}")
    print(f"Number of possible 5-pizza sets: {count_5_pizza_sets}")
    print("---")
    print(f"Total possible sets: {count_3_pizza_sets} + {count_4_pizza_sets} + {count_5_pizza_sets} = {total_sets}")
    
    return total_sets

# Run the calculation and store the final answer.
final_answer = solve_pizza_combinations()