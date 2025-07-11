import math

def find_pizza_combinations():
    """
    Calculates the number of possible sets of 3, 4, and 5 pizza sizes
    based on a specific set of constraints.
    """
    
    # --- Step 1 & 2: Find and Group Valid Pizzas ---
    
    # Dictionary to store valid diameters for each possible rounded slice area.
    # Key: rounded slice area (float), Value: set of diameters (int).
    diameters_by_area = {}

    # Constraint values
    MIN_DIAMETER = 1
    MAX_DIAMETER = 20
    MIN_SLICE_AREA = 5.0
    MAX_SLICE_AREA = 20.0

    # Iterate through all possible diameters (1 to 20).
    for diameter in range(MIN_DIAMETER, MAX_DIAMETER + 1):
        # Iterate through possible numbers of slices (multiples of 4).
        # A reasonable upper bound for slices is determined by the max pizza area
        # divided by the min slice area: (pi*10^2)/5 ~= 63. We check up to 100.
        for num_slices in range(4, 101, 4):
            pizza_total_area = math.pi * (diameter / 2) ** 2
            slice_area = pizza_total_area / num_slices
            rounded_area = round(slice_area, 1)

            # Check if this configuration meets the slice area constraint.
            if MIN_SLICE_AREA <= rounded_area <= MAX_SLICE_AREA:
                # If this is the first time we see this area, create a new set.
                if rounded_area not in diameters_by_area:
                    diameters_by_area[rounded_area] = set()
                # Add the diameter to the set for this slice area.
                diameters_by_area[rounded_area].add(diameter)

    # --- Step 3: Calculate Combinations ---

    count_sets_of_3 = 0
    count_sets_of_4 = 0
    count_sets_of_5 = 0

    # For each group of diameters, calculate the number of combinations.
    for area, diameters in diameters_by_area.items():
        n = len(diameters)
        
        # Combinations for sets of 3
        if n >= 3:
            count_sets_of_3 += math.comb(n, 3)
            
        # Combinations for sets of 4
        if n >= 4:
            count_sets_of_4 += math.comb(n, 4)
            
        # Combinations for sets of 5
        if n >= 5:
            count_sets_of_5 += math.comb(n, 5)

    # --- Step 4: Sum Totals and Print Results ---
    
    total_sets = count_sets_of_3 + count_sets_of_4 + count_sets_of_5

    print("To find the total number of possible pizza size sets, we first calculate the number of valid combinations for each set size and then sum them up.")
    print("\nCalculation of possible sets:")
    print(f"Sets of three pizza sizes: {count_sets_of_3}")
    print(f"Sets of four pizza sizes:  {count_sets_of_4}")
    print(f"Sets of five pizza sizes:  {count_sets_of_5}")

    print("\nThe final equation for the total number of sets is:")
    print(f"{count_sets_of_3} + {count_sets_of_4} + {count_sets_of_5} = {total_sets}")
    
    print("\nTherefore, the total number of possible sets is:")
    print(total_sets)

# Run the calculation.
find_pizza_combinations()
<<<14143>>>