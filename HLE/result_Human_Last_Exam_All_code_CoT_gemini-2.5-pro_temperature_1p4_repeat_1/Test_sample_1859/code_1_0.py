import math
from collections import defaultdict

def solve_pizza_combinations():
    """
    Calculates the number of possible sets of 3, 4, and 5 pizza sizes
    that meet the specified franchising constraints.
    """
    # --- Constraints ---
    MAX_DIAMETER = 20
    MIN_SLICE_AREA = 5.0
    MAX_SLICE_AREA = 20.0
    # A safe upper limit for the number of slices per pizza
    MAX_SLICES_TO_CHECK = 100 

    # This dictionary will store sets of valid diameters for each rounded slice area.
    # e.g., {10.5: {10, 12, 14}, 12.1: {11, 13, 15, 17}}
    valid_pizzas_by_area = defaultdict(set)

    # Step 1 & 2: Find all valid pizzas and group them by slice area.
    # Iterate through all possible integer diameters from 1 to 20 inches.
    for diameter in range(1, MAX_DIAMETER + 1):
        # Iterate through possible slice counts (must be divisible by 4).
        for num_slices in range(4, MAX_SLICES_TO_CHECK + 1, 4):
            # Calculate the area of the entire pizza.
            total_area = math.pi * (diameter / 2)**2
            # Calculate the area of a single slice.
            slice_area = total_area / num_slices
            # Round the slice area to one decimal place as required.
            rounded_slice_area = round(slice_area, 1)

            # Check if this pizza configuration meets the slice area constraint.
            if MIN_SLICE_AREA <= rounded_slice_area <= MAX_SLICE_AREA:
                # If it's a valid pizza, add its diameter to the set for its
                # corresponding slice area. Using a set prevents duplicate diameters.
                valid_pizzas_by_area[rounded_slice_area].add(diameter)

    # Step 3: For each group, count the combinations of 3, 4, and 5.
    total_combinations_of_3 = 0
    total_combinations_of_4 = 0
    total_combinations_of_5 = 0

    # Iterate through each set of diameters found for a common slice area.
    for diameters_set in valid_pizzas_by_area.values():
        n = len(diameters_set)
        
        # Calculate combinations for 3-pizza sets, if possible.
        if n >= 3:
            total_combinations_of_3 += math.comb(n, 3)
        # Calculate combinations for 4-pizza sets, if possible.
        if n >= 4:
            total_combinations_of_4 += math.comb(n, 4)
        # Calculate combinations for 5-pizza sets, if possible.
        if n >= 5:
            total_combinations_of_5 += math.comb(n, 5)

    # Step 4: Sum the results and print the breakdown.
    final_total = total_combinations_of_3 + total_combinations_of_4 + total_combinations_of_5

    print(f"Number of 3-pizza sets: {total_combinations_of_3}")
    print(f"Number of 4-pizza sets: {total_combinations_of_4}")
    print(f"Number of 5-pizza sets: {total_combinations_of_5}")
    print("---")
    print("The final equation for the total number of possible sets is:")
    print(f"{total_combinations_of_3} (3-pizza sets) + {total_combinations_of_4} (4-pizza sets) + {total_combinations_of_5} (5-pizza sets) = {final_total}")
    print(f"\nTotal possible franchising sets: {final_total}")


solve_pizza_combinations()