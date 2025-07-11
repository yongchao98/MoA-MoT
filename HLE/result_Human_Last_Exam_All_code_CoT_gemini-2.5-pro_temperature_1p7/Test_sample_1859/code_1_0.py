import math
from collections import defaultdict

def solve_pizza_combinations():
    """
    Calculates the number of possible sets of 3, 4, and 5 pizza sizes
    that meet the specified franchising criteria.
    """
    # --- Constraints ---
    DIAMETER_MAX = 20
    MIN_SLICE_AREA = 5.0
    MAX_SLICE_AREA = 20.0
    # A reasonable upper bound for slice count. Max pizza area / min slice area
    # is ~314 / 5 = ~63, so 100 is a safe upper limit.
    MAX_SLICES = 100

    # Dictionary to store diameters grouped by their valid rounded slice area
    # Key: rounded_slice_area, Value: list of diameters
    slice_area_to_diameters = defaultdict(list)

    # --- Step 1 & 2: Identify valid pizzas and group them by slice area ---
    # Iterate through all possible integer diameters
    for diameter in range(1, DIAMETER_MAX + 1):
        radius = diameter / 2.0
        total_area = math.pi * (radius**2)
        # Iterate through possible slice counts (must be divisible by 4)
        for slices in range(4, MAX_SLICES + 1, 4):
            area_per_slice = total_area / slices
            rounded_area = round(area_per_slice, 1)

            # Check if the rounded slice area is within the allowed range
            if MIN_SLICE_AREA <= rounded_area <= MAX_SLICE_AREA:
                slice_area_to_diameters[rounded_area].append(diameter)

    # --- Step 3: Count combinations for sets of 3, 4, and 5 ---
    count_sets_of_3 = 0
    count_sets_of_4 = 0
    count_sets_of_5 = 0

    # Iterate through each slice area and its corresponding list of possible diameters
    for area, diameters in slice_area_to_diameters.items():
        # A set must consist of distinct pizza sizes (diameters)
        unique_diameters = sorted(list(set(diameters)))
        n = len(unique_diameters)

        # Calculate combinations for sets of 3, if there are at least 3 unique diameters
        if n >= 3:
            count_sets_of_3 += math.comb(n, 3)

        # Calculate combinations for sets of 4, if there are at least 4 unique diameters
        if n >= 4:
            count_sets_of_4 += math.comb(n, 4)

        # Calculate combinations for sets of 5, if there are at least 5 unique diameters
        if n >= 5:
            count_sets_of_5 += math.comb(n, 5)

    # --- Step 4: Sum the results and print the final output ---
    total_sets = count_sets_of_3 + count_sets_of_4 + count_sets_of_5

    print(f"Number of possible sets with 3 pizza sizes: {count_sets_of_3}")
    print(f"Number of possible sets with 4 pizza sizes: {count_sets_of_4}")
    print(f"Number of possible sets with 5 pizza sizes: {count_sets_of_5}")
    print(f"Total number of possible sets is {count_sets_of_3} + {count_sets_of_4} + {count_sets_of_5} = {total_sets}")

solve_pizza_combinations()
<<<2836>>>