import math
from collections import defaultdict

def calculate_pizza_combinations():
    """
    Calculates the number of possible sets of 3, 4, and 5 pizza sizes
    based on a set of constraints.
    """
    
    # Step 1 & 2: Find all valid (diameter, slice_area) pairs and group by slice_area.
    slice_area_to_diameters = defaultdict(set)
    min_diameter = 1
    max_diameter = 20

    for diameter in range(min_diameter, max_diameter + 1):
        total_area = math.pi * (diameter / 2) ** 2
        
        # Determine the valid range for the number of slices (n)
        # 5 <= (total_area / n) <= 20  =>  total_area / 20 <= n <= total_area / 5
        # Also, n must be a multiple of 4.
        if total_area < 5 * 4: # Minimum area for 4 slices of 5 sq. in.
            continue
            
        min_n = math.ceil(total_area / 20)
        max_n = math.floor(total_area / 5)
        
        # Start checking from the first multiple of 4 within the valid range.
        start_n = ((int(min_n) + 3) // 4) * 4
        
        for num_slices in range(start_n, int(max_n) + 1, 4):
            if num_slices == 0:
                continue

            slice_area = total_area / num_slices
            rounded_slice_area = round(slice_area, 1)

            if 5.0 <= rounded_slice_area <= 20.0:
                slice_area_to_diameters[rounded_slice_area].add(diameter)

    # Step 3, 4, 5: Calculate combinations and sum the results.
    parts_3 = []
    parts_4 = []
    parts_5 = []

    # Iterate through sorted slice areas for deterministic output
    for area in sorted(slice_area_to_diameters.keys()):
        diameters = slice_area_to_diameters[area]
        k = len(diameters)

        if k >= 3:
            comb_3 = math.comb(k, 3)
            parts_3.append(str(comb_3))
        
        if k >= 4:
            comb_4 = math.comb(k, 4)
            parts_4.append(str(comb_4))
            
        if k >= 5:
            comb_5 = math.comb(k, 5)
            parts_5.append(str(comb_5))

    # Calculate and print the final results
    total_combinations_3 = sum(int(p) for p in parts_3)
    print("Combinations for sets of 3 pizzas:")
    print(f"{' + '.join(parts_3)} = {total_combinations_3}")

    total_combinations_4 = sum(int(p) for p in parts_4)
    print("\nCombinations for sets of 4 pizzas:")
    print(f"{' + '.join(parts_4)} = {total_combinations_4}")
    
    total_combinations_5 = sum(int(p) for p in parts_5)
    print("\nCombinations for sets of 5 pizzas:")
    print(f"{' + '.join(parts_5)} = {total_combinations_5}")
    
    grand_total = total_combinations_3 + total_combinations_4 + total_combinations_5
    print("\n--------------------------------------------------")
    print(f"Total possible sets of three, four, and five pizza sizes: {grand_total}")
    print("--------------------------------------------------")

if __name__ == '__main__':
    calculate_pizza_combinations()