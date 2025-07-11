def solve_puzzle():
    """
    Solves the matrix puzzle by calculating the missing elements and their sum.
    """

    # --- Step 1: Calculate the middle triplet of Row 3 ---
    # The triplet to the left is [7, 2, 9]
    x1, y1, z1 = 7, 2, 9
    
    # Check the condition x + y
    if (x1 + y1) > 10:
        # This block is not executed for this step
        next_x1 = (x1 * 3 - y1) % 12
        next_y1 = (y1 * 2 + 4) % 12
        next_z1 = (z1 + x1) % 12
    else: # x + y <= 10
        next_x1 = (x1 * 2 + y1) % 12
        next_y1 = (y1 * 3 - 2) % 12
        next_z1 = (z1 * 2) % 12
        
    missing_triplet_1 = [next_x1, next_y1, next_z1]
    
    # --- Step 2: Calculate the right triplet of Row 3 ---
    # The triplet to the left is the newly found middle triplet
    x2, y2, z2 = missing_triplet_1[0], missing_triplet_1[1], missing_triplet_1[2]

    # Check the condition x + y
    if (x2 + y2) > 10:
        # This block is not executed for this step
        next_x2 = (x2 * 3 - y2) % 12
        next_y2 = (y2 * 2 + 4) % 12
        next_z2 = (z2 + x2) % 12
    else: # x + y <= 10
        next_x2 = (x2 * 2 + y2) % 12
        next_y2 = (y2 * 3 - 2) % 12
        next_z2 = (z2 * 2) % 12
        
    missing_triplet_2 = [next_x2, next_y2, next_z2]
    
    # --- Step 3: Sum the missing elements and print the result ---
    all_missing_elements = missing_triplet_1 + missing_triplet_2
    total_sum = sum(all_missing_elements)
    
    # Create the equation string
    equation_str = " + ".join(map(str, all_missing_elements))
    
    print("The missing triplets are:")
    print(f"Middle triplet of Row 3: {missing_triplet_1}")
    print(f"Right triplet of Row 3: {missing_triplet_2}")
    print("\nThe sum of the missing elements is calculated as follows:")
    print(f"{equation_str} = {total_sum}")

solve_puzzle()
<<<24>>>