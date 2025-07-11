import math
from collections import defaultdict

def solve_pizza_combinations():
    """
    Calculates the number of possible sets of 3, 4, and 5 pizza sizes
    that meet the given constraints.
    """
    # This dictionary will map a rounded slice area to a set of possible diameters.
    pizzas_by_area = defaultdict(set)

    # Define constraints
    MAX_DIAMETER = 20
    MIN_SLICE_AREA = 5.0
    MAX_SLICE_AREA = 20.0
    
    # Step 1: Generate all possible valid individual pizzas and group them by slice area.
    
    # Iterate through all possible integer diameters from 1 to 20 inches.
    for diameter in range(1, MAX_DIAMETER + 1):
        pizza_area = math.pi * (diameter / 2)**2
        
        # Iterate through possible numbers of slices (must be divisible by 4).
        # The maximum possible number of slices is for the largest pizza (20")
        # and the smallest slice area (5.0 sq in), which is approx. 62.8.
        # So, we check multiples of 4 up to 60.
        for num_slices in range(4, 61, 4):
            # Calculate the area of a single slice
            slice_area = pizza_area / num_slices
            
            # Check if the actual slice area is within the allowed range
            if MIN_SLICE_AREA <= slice_area <= MAX_SLICE_AREA:
                # Round the slice area to one decimal place as per the rule
                rounded_area = round(slice_area, 1)
                
                # Add the diameter to the set for this rounded area.
                # Using a set automatically handles duplicates.
                pizzas_by_area[rounded_area].add(diameter)

    # Step 2: Calculate combinations for each slice area.
    
    total_3_pizza_sets = 0
    total_4_pizza_sets = 0
    total_5_pizza_sets = 0

    # Iterate through the groups of pizzas categorized by their rounded slice area
    for rounded_area, diameters in pizzas_by_area.items():
        num_diameters = len(diameters)
        
        # Calculate combinations for sets of 3 pizza sizes
        if num_diameters >= 3:
            total_3_pizza_sets += math.comb(num_diameters, 3)
            
        # Calculate combinations for sets of 4 pizza sizes
        if num_diameters >= 4:
            total_4_pizza_sets += math.comb(num_diameters, 4)
            
        # Calculate combinations for sets of 5 pizza sizes
        if num_diameters >= 5:
            total_5_pizza_sets += math.comb(num_diameters, 5)

    # Step 3: Sum the totals and print the final result.
    
    total_sets = total_3_pizza_sets + total_4_pizza_sets + total_5_pizza_sets

    print(f"Possible sets of 3 pizza sizes: {total_3_pizza_sets}")
    print(f"Possible sets of 4 pizza sizes: {total_4_pizza_sets}")
    print(f"Possible sets of 5 pizza sizes: {total_5_pizza_sets}")
    print("\nTotal possible sets:")
    print(f"{total_3_pizza_sets} + {total_4_pizza_sets} + {total_5_pizza_sets} = {total_sets}")
    
    return total_sets

# Execute the function to find the answer
total_combinations = solve_pizza_combinations()
print(f"\n<<<{total_combinations}>>>")