import math

def solve_pizza_combinations():
    """
    Calculates the number of possible sets of 3, 4, and 5 pizza sizes
    that meet the given franchising constraints.
    """
    
    # --- Step 1 & 2: Generate and filter for valid pizzas, then group them ---
    # A dictionary to group unique diameters by their resulting slice area.
    # The key is the rounded_slice_area, and the value is a set of diameters.
    pizzas_by_area = {}

    # Iterate through all possible integer diameters from 1 to 20 inches.
    for diameter in range(1, 21):
        # Iterate through possible numbers of slices (must be a multiple of 4).
        # A maximum of 100 slices is a safe upper bound, as the max slice area
        # constraint limits the required number of slices.
        for slices in range(4, 101, 4):
            # Calculate the area of a single slice.
            pizza_area = math.pi * (diameter / 2)**2
            slice_area = pizza_area / slices
            
            # Round the slice area to one decimal place, as per requirements.
            rounded_slice_area = round(slice_area, 1)

            # Check if the rounded slice area falls within the allowed range (5.0 to 20.0).
            if 5.0 <= rounded_slice_area <= 20.0:
                # If this slice area hasn't been seen before, initialize it with an empty set.
                if rounded_slice_area not in pizzas_by_area:
                    pizzas_by_area[rounded_slice_area] = set()
                # Add the current diameter to the set for this slice area.
                # Using a set automatically handles duplicates.
                pizzas_by_area[rounded_slice_area].add(diameter)

    # --- Step 3: Calculate the number of combinations for each set size ---
    
    total_sets_of_3 = 0
    total_sets_of_4 = 0
    total_sets_of_5 = 0

    # Iterate through the unique slice areas we found.
    for area in pizzas_by_area:
        # Get the number of unique diameters that can produce this slice area.
        num_diameters = len(pizzas_by_area[area])

        # Calculate how many sets of 3, 4, and 5 can be formed from these diameters.
        # math.comb(n, k) calculates "n choose k".
        
        # Add combinations of 3 pizzas if there are at least 3 diameters available.
        if num_diameters >= 3:
            total_sets_of_3 += math.comb(num_diameters, 3)
            
        # Add combinations of 4 pizzas if there are at least 4 diameters available.
        if num_diameters >= 4:
            total_sets_of_4 += math.comb(num_diameters, 4)
            
        # Add combinations of 5 pizzas if there are at least 5 diameters available.
        if num_diameters >= 5:
            total_sets_of_5 += math.comb(num_diameters, 5)
            
    # --- Step 4: Print the results ---
    
    # Calculate the total number of sets.
    total_combinations = total_sets_of_3 + total_sets_of_4 + total_sets_of_5

    print(f"Possible sets of 3 pizza sizes: {total_sets_of_3}")
    print(f"Possible sets of 4 pizza sizes: {total_sets_of_4}")
    print(f"Possible sets of 5 pizza sizes: {total_sets_of_5}")
    print("---")
    print(f"Final Equation: {total_sets_of_3} + {total_sets_of_4} + {total_sets_of_5} = {total_combinations}")
    print(f"Total possible sets: {total_combinations}")

# Execute the main function to find and print the solution.
solve_pizza_combinations()
<<<1208>>>