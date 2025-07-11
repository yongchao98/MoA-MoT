import math
from collections import defaultdict

def find_pizza_combinations():
    """
    Calculates the number of possible sets of 3, 4, and 5 pizza sizes
    based on a given set of constraints.
    """

    # Step 1 & 2: Identify valid pizzas and group them by slice area.
    # The dictionary will map a slice area to a set of diameters that can produce it.
    area_to_diameters = defaultdict(set)
    
    # Constraints
    max_diameter = 20
    min_slice_area = 5.0
    max_slice_area = 20.0
    
    # Iterate through all possible diameters (1 to 20 inches).
    for d in range(1, max_diameter + 1):
        # Iterate through possible numbers of slices (must be divisible by 4).
        # A slice count above 60 is not feasible as the slice area would be too small.
        for n in range(4, 61, 4):
            # Calculate the area of a single slice.
            # Pizza Area = pi * r^2 = pi * (d/2)^2 = (pi * d^2) / 4
            # Slice Area = Pizza Area / n
            slice_area = (math.pi * d**2) / (4 * n)
            
            # Round the area to one decimal place as required.
            rounded_area = round(slice_area, 1)
            
            # Check if the slice area is within the allowed range.
            if min_slice_area <= rounded_area <= max_slice_area:
                # If valid, add this diameter to the set for this specific slice area.
                area_to_diameters[rounded_area].add(d)

    def combinations(n, k):
        """
        Calculates the number of combinations (n choose k).
        Returns 0 if k > n.
        """
        if k < 0 or k > n:
            return 0
        # Use math.comb for efficiency and simplicity.
        return math.comb(n, k)

    # Step 3: Calculate the number of combinations for each group.
    count_3_pizza_sets = 0
    count_4_pizza_sets = 0
    count_5_pizza_sets = 0

    # Iterate through the dictionary of grouped diameters.
    for area in area_to_diameters:
        num_available_diameters = len(area_to_diameters[area])
        
        # Calculate how many sets of 3, 4, and 5 can be made.
        count_3_pizza_sets += combinations(num_available_diameters, 3)
        count_4_pizza_sets += combinations(num_available_diameters, 4)
        count_5_pizza_sets += combinations(num_available_diameters, 5)

    # Step 4: Present the final results.
    print(f"Possible sets of three pizza sizes: {count_3_pizza_sets}")
    print(f"Possible sets of four pizza sizes: {count_4_pizza_sets}")
    print(f"Possible sets of five pizza sizes: {count_5_pizza_sets}")
    
    total_sets = count_3_pizza_sets + count_4_pizza_sets + count_5_pizza_sets
    print(f"\nTotal possible sets: {count_3_pizza_sets} + {count_4_pizza_sets} + {count_5_pizza_sets} = {total_sets}")

# Execute the function to find and print the answer.
find_pizza_combinations()
<<<2585>>>