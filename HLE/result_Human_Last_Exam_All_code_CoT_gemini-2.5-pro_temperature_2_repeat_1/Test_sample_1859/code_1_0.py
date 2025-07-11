import math
from collections import defaultdict

def find_pizza_combinations():
    """
    Calculates the number of possible sets of 3, 4, and 5 pizza sizes
    that meet the specified franchising criteria.
    """
    # Define problem constraints
    MAX_DIAMETER = 20
    MIN_SLICE_AREA = 5.0
    MAX_SLICE_AREA = 20.0
    # Heuristic maximum for slices; a 20" pizza with 5 sq.in. slices has ~63 slices.
    MAX_SLICES_TO_CHECK = 80

    # A dictionary to hold sets of valid diameters for each possible slice area.
    # e.g., {9.8: {7, 10, 13, 15, 17}, 14.1: {8, 12, 16}}
    diameters_by_area = defaultdict(set)

    # Step 1: Find all valid pizza diameters for each possible slice area.
    for diameter in range(1, MAX_DIAMETER + 1):
        for n_slices in range(4, MAX_SLICES_TO_CHECK + 1, 4):
            pizza_area = math.pi * (diameter / 2) ** 2
            current_slice_area = pizza_area / n_slices

            # The area of each slice must be rounded to one decimal place.
            rounded_slice_area = round(current_slice_area, 1)

            # Check if the calculated slice area falls within the allowed range.
            if MIN_SLICE_AREA <= rounded_slice_area <= MAX_SLICE_AREA:
                diameters_by_area[rounded_slice_area].add(diameter)

    # Step 2: Calculate the number of combinations for sets of 3, 4, and 5.
    total_sets_3 = 0
    total_sets_4 = 0
    total_sets_5 = 0

    for area in diameters_by_area:
        num_valid_diameters = len(diameters_by_area[area])

        # For each common slice area, calculate how many sets can be formed.
        if num_valid_diameters >= 3:
            total_sets_3 += math.comb(num_valid_diameters, 3)
        if num_valid_diameters >= 4:
            total_sets_4 += math.comb(num_valid_diameters, 4)
        if num_valid_diameters >= 5:
            total_sets_5 += math.comb(num_valid_diameters, 5)

    # Step 3: Sum the results and print the final calculation.
    total_combinations = total_sets_3 + total_sets_4 + total_sets_5

    print(f"Number of possible sets with 3 pizza sizes: {total_sets_3}")
    print(f"Number of possible sets with 4 pizza sizes: {total_sets_4}")
    print(f"Number of possible sets with 5 pizza sizes: {total_sets_5}")
    print("\n---")
    print("Total possible sets calculation:")
    print(f"{total_sets_3} (3-pizza sets) + {total_sets_4} (4-pizza sets) + {total_sets_5} (5-pizza sets) = {total_combinations}")


if __name__ == "__main__":
    find_pizza_combinations()
<<<1136>>>