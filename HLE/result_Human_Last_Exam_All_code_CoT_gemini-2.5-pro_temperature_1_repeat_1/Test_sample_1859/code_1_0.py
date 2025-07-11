import math

def find_pizza_combinations():
    """
    Calculates the number of possible sets of 3, 4, and 5 pizza sizes
    based on a given set of constraints.
    """
    # This dictionary will map a rounded slice area to a set of valid diameters.
    # Using a set for diameters automatically handles duplicates.
    valid_diameters_by_area = {}

    max_diameter = 20
    # The number of slices (N) can be estimated: N ≈ (π * d²) / (4 * S).
    # The maximum N occurs at max diameter (20) and min slice area (5),
    # which is N ≈ (π * 20²) / (4 * 5) ≈ 62.8.
    # A loop up to 100 slices is more than sufficient.
    max_slices = 100

    # Step 1: Iterate through all possible pizza configurations to find valid diameters for each possible slice area.
    for diameter in range(1, max_diameter + 1):
        # The number of slices must be a multiple of 4.
        for slices in range(4, max_slices + 1, 4):
            # Calculate the precise area of a single slice.
            # Pizza Area = π * r² = π * (diameter/2)²
            slice_area = (math.pi * diameter**2) / (4 * slices)
            
            # Round the slice area to one decimal place as required.
            rounded_area = round(slice_area, 1)
            
            # Check if the rounded area is within the allowed range [5.0, 20.0].
            if 5.0 <= rounded_area <= 20.0:
                # If this area is new, initialize a set to store its diameters.
                if rounded_area not in valid_diameters_by_area:
                    valid_diameters_by_area[rounded_area] = set()
                # Add the current diameter to the set for this slice area.
                valid_diameters_by_area[rounded_area].add(diameter)

    # Step 2: Calculate the number of combinations for each group of diameters.
    total_combinations_3 = 0
    total_combinations_4 = 0
    total_combinations_5 = 0

    # Iterate through the dictionary of {area: {diameters}}.
    for area, diameters_set in valid_diameters_by_area.items():
        num_possible_diameters = len(diameters_set)
        
        # Calculate combinations of 3 pizza sizes if there are at least 3 diameters.
        if num_possible_diameters >= 3:
            total_combinations_3 += math.comb(num_possible_diameters, 3)
            
        # Calculate combinations of 4 pizza sizes if there are at least 4 diameters.
        if num_possible_diameters >= 4:
            total_combinations_4 += math.comb(num_possible_diameters, 4)

        # Calculate combinations of 5 pizza sizes if there are at least 5 diameters.
        if num_possible_diameters >= 5:
            total_combinations_5 += math.comb(num_possible_diameters, 5)
            
    # Step 3: Sum the results and print the final breakdown.
    total_combinations = total_combinations_3 + total_combinations_4 + total_combinations_5

    print(f"Number of 3-pizza-size combinations: {total_combinations_3}")
    print(f"Number of 4-pizza-size combinations: {total_combinations_4}")
    print(f"Number of 5-pizza-size combinations: {total_combinations_5}")
    print(f"Total possible sets of pizza sizes = {total_combinations_3} + {total_combinations_4} + {total_combinations_5} = {total_combinations}")

if __name__ == '__main__':
    find_pizza_combinations()