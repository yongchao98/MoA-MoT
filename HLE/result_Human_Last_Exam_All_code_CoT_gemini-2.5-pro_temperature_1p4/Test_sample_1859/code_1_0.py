import math

def solve_pizza_combinations():
    """
    Calculates the number of possible sets of three, four, and five pizza sizes
    that meet the given franchising constraints.
    """

    # This dictionary will store the sets of valid diameters for each possible rounded slice area.
    # The key is the rounded_slice_area (float), and the value is a set of diameters (int).
    valid_pizzas_by_area = {}

    # Define the problem constraints
    MAX_DIAMETER = 20
    MIN_SLICE_AREA = 5.0
    MAX_SLICE_AREA = 20.0

    # Step 1 & 2: Iterate through all possible pizzas to find and filter valid ones.
    # A pizza is defined by its diameter and number of slices.
    
    # Iterate through all possible integer diameters from 1 to 20 inches.
    for diameter in range(1, MAX_DIAMETER + 1):
        # Iterate through a reasonable range of slice counts (must be multiples of 4).
        # A maximum of 100 slices is a safe upper bound for a 20-inch pizza.
        for num_slices in range(4, 101, 4):
            # Calculate the area of a single slice using the formula:
            # Slice Area = Pizza Area / Number of Slices
            # Pizza Area = pi * r^2 = pi * (diameter / 2)^2
            slice_area_raw = (math.pi * (diameter ** 2)) / (4 * num_slices)
            
            # Round the slice area to one decimal place as per the requirements.
            rounded_slice_area = round(slice_area_raw, 1)
            
            # Check if the calculated slice area is within the allowed range.
            if MIN_SLICE_AREA <= rounded_slice_area <= MAX_SLICE_AREA:
                # Step 3: Group the valid diameters by their common slice area.
                # If this slice area hasn't been seen before, initialize a new set for it.
                if rounded_slice_area not in valid_pizzas_by_area:
                    valid_pizzas_by_area[rounded_slice_area] = set()
                # Add the current diameter to the set for this specific slice area.
                valid_pizzas_by_area[rounded_slice_area].add(diameter)

    # Step 4: Calculate the number of combinations for sets of 3, 4, and 5.
    total_sets_of_3 = 0
    total_sets_of_4 = 0
    total_sets_of_5 = 0

    # Iterate through each slice area and its corresponding set of compatible diameters.
    for diameters_set in valid_pizzas_by_area.values():
        num_compatible_diameters = len(diameters_set)
        
        # Calculate combinations for sets of 3, if there are at least 3 compatible diameters.
        if num_compatible_diameters >= 3:
            total_sets_of_3 += math.comb(num_compatible_diameters, 3)
            
        # Calculate combinations for sets of 4, if there are at least 4 compatible diameters.
        if num_compatible_diameters >= 4:
            total_sets_of_4 += math.comb(num_compatible_diameters, 4)
            
        # Calculate combinations for sets of 5, if there are at least 5 compatible diameters.
        if num_compatible_diameters >= 5:
            total_sets_of_5 += math.comb(num_compatible_diameters, 5)

    # Step 5: Sum and report the final results.
    total_combinations = total_sets_of_3 + total_sets_of_4 + total_sets_of_5

    print(f"Possible sets of three pizza sizes: {total_sets_of_3}")
    print(f"Possible sets of four pizza sizes: {total_sets_of_4}")
    print(f"Possible sets of five pizza sizes: {total_sets_of_5}")
    print(f"\nTotal possible sets: {total_sets_of_3} + {total_sets_of_4} + {total_sets_of_5} = {total_combinations}")
    
    return total_combinations

# Run the calculation and store the final answer.
final_answer = solve_pizza_combinations()
# The final answer needs to be enclosed in "<<< >>>" as requested.
# However, the user wants to see the printed output from the function, 
# so the answer block will just contain the final numerical result.
# The code above already prints the detailed breakdown.
final_answer_str = f"<<<{final_answer}>>>"
print(final_answer_str)