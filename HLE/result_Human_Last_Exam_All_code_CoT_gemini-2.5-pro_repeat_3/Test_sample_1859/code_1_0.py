import math

def solve_pizza_combinations():
    """
    Calculates the number of possible sets of 3, 4, and 5 pizza sizes
    that meet the given constraints.
    """
    # Step 1 & 2: Find all valid pizza configurations and group them by slice area.
    # A dictionary is used to store a set of possible diameters for each valid slice area.
    pizzas_by_slice_area = {}
    min_diameter = 1
    max_diameter = 20
    min_slice_area = 5.0
    max_slice_area = 20.0

    # Iterate through all possible integer diameters from 1 to 20 inches.
    for diameter in range(min_diameter, max_diameter + 1):
        # Iterate through possible numbers of slices (multiples of 4).
        # A generous upper bound (100) is used to ensure all possibilities are covered.
        for num_slices in range(4, 101, 4):
            # Calculate the area of the pizza and then the area of a single slice.
            radius = diameter / 2
            pizza_area = math.pi * (radius ** 2)
            slice_area = pizza_area / num_slices
            
            # Round the slice area to one decimal place as required.
            rounded_slice_area = round(slice_area, 1)

            # Check if the rounded slice area is within the allowed range [5.0, 20.0].
            if min_slice_area <= rounded_slice_area <= max_slice_area:
                # If this slice area hasn't been seen, initialize a new set for it.
                if rounded_slice_area not in pizzas_by_slice_area:
                    pizzas_by_slice_area[rounded_slice_area] = set()
                # Add the current diameter to the set for this specific slice area.
                pizzas_by_slice_area[rounded_slice_area].add(diameter)

    # Step 3 & 4: Calculate combinations for sets of 3, 4, and 5, and sum the totals.
    combinations_of_3 = 0
    combinations_of_4 = 0
    combinations_of_5 = 0

    # Iterate through the groups of diameters for each common slice area.
    for slice_area, diameters in pizzas_by_slice_area.items():
        num_possible_diameters = len(diameters)

        # Calculate combinations of 3 if there are at least 3 diameters available.
        if num_possible_diameters >= 3:
            combinations_of_3 += math.comb(num_possible_diameters, 3)
        
        # Calculate combinations of 4 if there are at least 4 diameters available.
        if num_possible_diameters >= 4:
            combinations_of_4 += math.comb(num_possible_diameters, 4)

        # Calculate combinations of 5 if there are at least 5 diameters available.
        if num_possible_diameters >= 5:
            combinations_of_5 += math.comb(num_possible_diameters, 5)

    # Print the final counts for each set size.
    print(f"Number of possible sets with 3 pizza sizes: {combinations_of_3}")
    print(f"Number of possible sets with 4 pizza sizes: {combinations_of_4}")
    print(f"Number of possible sets with 5 pizza sizes: {combinations_of_5}")
    
    total_combinations = combinations_of_3 + combinations_of_4 + combinations_of_5
    print(f"Total number of possible sets: {total_combinations}")


# Execute the main function to find and print the solution.
solve_pizza_combinations()
<<<3645>>>