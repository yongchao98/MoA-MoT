import math
from collections import defaultdict

def find_pizza_combinations():
    """
    Calculates the number of possible sets of 3, 4, and 5 pizza sizes
    that meet a series of franchising constraints.

    The constraints are:
    - Diameters are integers up to 20 inches.
    - The number of slices is divisible by 4.
    - Slice area is consistent across a set (rounded to one decimal).
    - Slice area is between 5.0 and 20.0 sq. inches, inclusive.
    """
    
    # Define constraints
    MAX_DIAMETER = 20
    MIN_SLICE_AREA = 5.0
    MAX_SLICE_AREA = 20.0
    
    # This dictionary will store lists of diameters for each valid slice area.
    # e.g., {15.7: [10, 14, 18], 12.6: [8, 12, 16]}
    area_to_diameters = defaultdict(list)
    
    # Step 1 & 2: Find all valid individual pizzas and group them by slice area.
    
    # Iterate through all possible integer diameters from 1 to 20 inches.
    for diameter in range(1, MAX_DIAMETER + 1):
        pizza_area = math.pi * (diameter / 2)**2
        
        # Iterate through a reasonable range of slice counts (must be divisible by 4).
        # Maximum possible slices for a 20-inch pizza with 5 sq. in. slices is ~62.
        # We can safely check up to 64 slices.
        for num_slices in range(4, 65, 4):
            # Calculate slice area and round to one decimal place per the rule.
            slice_area = round(pizza_area / num_slices, 1)
            
            # Check if this configuration is valid.
            if MIN_SLICE_AREA <= slice_area <= MAX_SLICE_AREA:
                # Add the diameter to the list for this specific slice area.
                # We prevent duplicates in case a diameter with a different
                # slice count produces the same rounded area.
                if diameter not in area_to_diameters[slice_area]:
                    area_to_diameters[slice_area].append(diameter)

    # Step 3: Calculate the number of combinations for sets of 3, 4, and 5.
    
    combinations_of_3 = 0
    combinations_of_4 = 0
    combinations_of_5 = 0
    
    # Iterate through each group of compatible diameters.
    for diameters in area_to_diameters.values():
        n = len(diameters)
        
        # Calculate how many sets of 3 can be made if we have at least 3 options.
        if n >= 3:
            combinations_of_3 += math.comb(n, 3)
            
        # Calculate how many sets of 4 can be made if we have at least 4 options.
        if n >= 4:
            combinations_of_4 += math.comb(n, 4)
            
        # Calculate how many sets of 5 can be made if we have at least 5 options.
        if n >= 5:
            combinations_of_5 += math.comb(n, 5)
            
    # Step 4: Sum the results and display the final equation.
    
    total_sets = combinations_of_3 + combinations_of_4 + combinations_of_5
    
    print(f"Possible sets with 3 pizza sizes: {combinations_of_3}")
    print(f"Possible sets with 4 pizza sizes: {combinations_of_4}")
    print(f"Possible sets with 5 pizza sizes: {combinations_of_5}")
    print("-" * 40)
    print(f"Total possible sets = {combinations_of_3} + {combinations_of_4} + {combinations_of_5} = {total_sets}")


if __name__ == "__main__":
    find_pizza_combinations()

<<<2242>>>