import math

def solve_pizza_combinations():
    """
    Calculates the number of possible sets of three, four, and five pizza sizes
    that meet a given set of constraints.
    """

    # Dictionary to store sets of diameters for each valid slice area.
    # Key: slice_area (float, rounded to 1 decimal place)
    # Value: set of diameters (int) that can produce this slice area.
    slice_area_to_diameters = {}

    # Step 1 & 2: Iterate through all possible pizza configurations to find valid ones.
    for diameter in range(1, 21):
        pizza_area = (math.pi / 4) * (diameter**2)
        n_slices = 4
        
        # We test slice counts (multiples of 4) until the slice area becomes too small.
        while True:
            # Calculate the raw and rounded slice area.
            raw_slice_area = pizza_area / n_slices
            rounded_area = round(raw_slice_area, 1)

            # If the rounded area drops below 5.0, any larger number of slices
            # will also be too small, so we can stop checking for this diameter.
            if rounded_area < 5.0:
                break
            
            # Step 3: If the rounded area is within the valid range [5.0, 20.0],
            # add the diameter to the set for that area.
            if rounded_area <= 20.0:
                if rounded_area not in slice_area_to_diameters:
                    slice_area_to_diameters[rounded_area] = set()
                slice_area_to_diameters[rounded_area].add(diameter)
            
            n_slices += 4

    # Step 4: Calculate the number of combinations for sets of 3, 4, and 5 pizzas.
    sets_of_3 = 0
    sets_of_4 = 0
    sets_of_5 = 0

    for area in slice_area_to_diameters:
        num_available_diameters = len(slice_area_to_diameters[area])
        
        # The math.comb(n, k) function automatically returns 0 if k > n.
        sets_of_3 += math.comb(num_available_diameters, 3)
        sets_of_4 += math.comb(num_available_diameters, 4)
        sets_of_5 += math.comb(num_available_diameters, 5)

    # Step 5: Sum the totals and print the final result.
    total_sets = sets_of_3 + sets_of_4 + sets_of_5

    print(f"Possible sets of 3 pizza sizes: {sets_of_3}")
    print(f"Possible sets of 4 pizza sizes: {sets_of_4}")
    print(f"Possible sets of 5 pizza sizes: {sets_of_5}")
    print(f"\nTotal possible sets: {sets_of_3} + {sets_of_4} + {sets_of_5} = {total_sets}")
    
# Run the function to solve the problem.
solve_pizza_combinations()