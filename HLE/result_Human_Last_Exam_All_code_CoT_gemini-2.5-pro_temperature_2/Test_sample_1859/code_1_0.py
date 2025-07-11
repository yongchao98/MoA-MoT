import math
from collections import defaultdict

def find_pizza_combinations():
    """
    Calculates the number of possible sets of three, four, and five pizza sizes 
    based on a series of restaurant-specific constraints.
    """
    # --- Step 1: Define Constraints ---
    MIN_DIAMETER = 1
    MAX_DIAMETER = 20
    MIN_SLICE_AREA = 5.0
    MAX_SLICE_AREA = 20.0
    # A safe upper bound for slices to check. For a 20-inch pizza, the maximum
    # number of 5 sq. inch slices is around 62, so 100 is plenty.
    MAX_SLICES_TO_CHECK = 100 

    # --- Step 2: Group valid diameters by common slice area ---
    # This dictionary will store a set of diameters for each slice area.
    # Example: {9.8: {10, 14, 16}, 11.3: {12, 17, 20}}
    valid_pizzas_by_area = defaultdict(set)

    # Iterate through each possible integer diameter.
    for diameter in range(MIN_DIAMETER, MAX_DIAMETER + 1):
        pizza_total_area = math.pi * (diameter / 2)**2
        
        # Iterate through possible numbers of slices (must be a multiple of 4).
        for num_slices in range(4, MAX_SLICES_TO_CHECK + 1, 4):
            # No need to check if the slice area is already too small.
            if pizza_total_area / num_slices < 4.95:
                break
            
            slice_area = pizza_total_area / num_slices
            
            # Round the slice area to one decimal place, as required.
            rounded_area = round(slice_area, 1)

            # Check if this pizza configuration meets the slice area constraint.
            if MIN_SLICE_AREA <= rounded_area <= MAX_SLICE_AREA:
                valid_pizzas_by_area[rounded_area].add(diameter)

    # --- Step 3: Calculate combinations for each set size ---
    total_sets_of_3 = 0
    total_sets_of_4 = 0
    total_sets_of_5 = 0

    # Iterate through each common slice area that was found.
    for area, diameters in valid_pizzas_by_area.items():
        num_valid_diameters = len(diameters)
        
        # Calculate how many sets of 3, 4, and 5 can be formed from the
        # available diameters for the current slice area.
        total_sets_of_3 += math.comb(num_valid_diameters, 3)
        total_sets_of_4 += math.comb(num_valid_diameters, 4)
        total_sets_of_5 += math.comb(num_valid_diameters, 5)

    # --- Step 4: Display the final results ---
    grand_total = total_sets_of_3 + total_sets_of_4 + total_sets_of_5

    print("To find the total number of possible pizza size sets, we calculate the combinations for sets of three, four, and five pizzas and sum them.\n")
    print(f"Possible sets with 3 pizza sizes: {total_sets_of_3}")
    print(f"Possible sets with 4 pizza sizes: {total_sets_of_4}")
    print(f"Possible sets with 5 pizza sizes: {total_sets_of_5}")
    print("\nThe final equation is:")
    print(f"{total_sets_of_3} + {total_sets_of_4} + {total_sets_of_5} = {grand_total}")
    print(f"\nIn total, there are {grand_total} possible sets you can create.")

if __name__ == '__main__':
    find_pizza_combinations()
    print("\n<<<1102>>>")