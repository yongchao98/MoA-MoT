import math

def solve_pizza_combinations():
    """
    Calculates the total number of possible sets of three, four, and five pizza sizes
    based on a given set of constraints.
    """
    # This dictionary will map a rounded slice area to a set of valid pizza diameters.
    # e.g., {9.8: {10, 14}, 12.6: {8, 11, 16}}
    area_to_diameters = {}

    # Constraints from the problem description
    MAX_DIAMETER = 20
    MIN_SLICE_AREA = 5.0
    MAX_SLICE_AREA = 20.0

    # Step 1 & 2: Generate all possible valid pizza configurations and filter by slice area.
    # Iterate through all possible integer diameters from 1 to 20.
    for d in range(1, MAX_DIAMETER + 1):
        # Iterate through possible numbers of slices (must be a multiple of 4).
        # A reasonable upper bound for slices is 64, as a larger number of slices with the
        # largest pizza would result in a slice area smaller than the minimum.
        for s in range(4, 68, 4):
            pizza_total_area = math.pi * (d / 2)**2
            slice_area = pizza_total_area / s

            # Round the area to one decimal place, as per requirements.
            rounded_area = round(slice_area, 1)

            # Check if the rounded area falls within the allowed range.
            if MIN_SLICE_AREA <= rounded_area <= MAX_SLICE_AREA:
                # Step 3: Group the valid diameters by their common rounded slice area.
                if rounded_area not in area_to_diameters:
                    area_to_diameters[rounded_area] = set()
                area_to_diameters[rounded_area].add(d)

    # Step 4: Calculate the number of combinations for each slice area group.
    sub_totals = []
    
    # Sort the keys (slice areas) for a consistent and readable output.
    sorted_areas = sorted(area_to_diameters.keys())

    for area in sorted_areas:
        diameters = area_to_diameters[area]
        n = len(diameters)
        
        current_area_combinations = 0
        
        # Calculate combinations for sets of 3, 4, or 5 pizza sizes.
        if n >= 3:
            current_area_combinations += math.comb(n, 3)
        if n >= 4:
            current_area_combinations += math.comb(n, 4)
        if n >= 5:
            current_area_combinations += math.comb(n, 5)

        if current_area_combinations > 0:
            sub_totals.append(current_area_combinations)

    # Step 5: Sum the sub-totals to get the final answer and print the results.
    total_possible_sets = sum(sub_totals)

    # To fulfill the requirement "output each number in the final equation"
    # we will show the sum of the subtotals from each group.
    equation_str = " + ".join(map(str, sub_totals))
    
    print("The total number of possible sets is calculated by summing the valid combinations from each group with a unique slice area.")
    print("\nFinal Calculation:")
    print(f"{equation_str} = {total_possible_sets}")
    
    print(f"\nIn total, there are {total_possible_sets} possible sets of three, four, and five pizza sizes that meet your criteria.")

solve_pizza_combinations()