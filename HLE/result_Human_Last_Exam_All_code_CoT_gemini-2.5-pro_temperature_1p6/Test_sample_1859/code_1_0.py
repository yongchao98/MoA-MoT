import math

def find_pizza_combinations():
    """
    Calculates the number of possible sets of 3, 4, and 5 pizza sizes
    that meet the given franchising constraints.
    """
    pizzas_by_area = {}
    max_diameter = 20
    min_slice_area = 5.0
    max_slice_area = 20.0
    
    # Step 1, 2 & 3: Generate, filter, and group all valid pizzas.
    # Iterate through all possible integer diameters from 1 to 20.
    for diameter in range(1, max_diameter + 1):
        # The number of slices must be a multiple of 4.
        # We can determine a reasonable upper limit for slices. For the largest
        # diameter (20) and smallest slice area (5), the max slices would be
        # (pi * 10^2) / 5 ~= 62.8. Iterating up to 64 is sufficient.
        for slices in range(4, 68, 4):
            # Calculate the precise area of a single slice.
            pizza_area = math.pi * (diameter / 2)**2
            slice_area = pizza_area / slices
            
            # Check if the actual slice area is within the valid range [5, 20].
            if min_slice_area <= slice_area <= max_slice_area:
                # Group pizzas by their slice area rounded to one decimal place.
                rounded_area = round(slice_area, 1)
                
                # Add the diameter to the set for this specific rounded area.
                if rounded_area not in pizzas_by_area:
                    pizzas_by_area[rounded_area] = set()
                pizzas_by_area[rounded_area].add(diameter)

    # Step 4: Calculate the number of combinations for each group.
    total_sets_of_3 = 0
    total_sets_of_4 = 0
    total_sets_of_5 = 0

    for area in pizzas_by_area:
        # 'k' is the number of different pizza diameters for a given slice area.
        k = len(pizzas_by_area[area])
        
        # Calculate combinations of 3, 4, and 5 from the k options.
        # math.comb(k, n) is 0 if n > k, so no need for checks.
        total_sets_of_3 += math.comb(k, 3)
        total_sets_of_4 += math.comb(k, 4)
        total_sets_of_5 += math.comb(k, 5)
            
    # Step 5: Sum the totals and print the final result.
    total_combinations = total_sets_of_3 + total_sets_of_4 + total_sets_of_5
    
    print("Based on the constraints, we can form the following sets:")
    print(f"Number of valid combinations for sets of 3 pizzas: {total_sets_of_3}")
    print(f"Number of valid combinations for sets of 4 pizzas: {total_sets_of_4}")
    print(f"Number of valid combinations for sets of 5 pizzas: {total_sets_of_5}")
    print("\nThe total number of possible sets is the sum of these combinations:")
    print(f"Total possible sets = {total_sets_of_3} + {total_sets_of_4} + {total_sets_of_5} = {total_combinations}")

# Execute the function to find the answer.
find_pizza_combinations()