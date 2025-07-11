def solve_puzzle():
    """
    This function calculates the missing elements in the matrix based on the provided rules
    and then calculates their sum.
    """

    # The first triplet in the third row is given.
    r3c1 = [7, 2, 9]
    print(f"Starting with the first known triplet in Row 3: {r3c1}")

    # --- Calculate the middle triplet of Row 3 (M[2][1]) ---
    x, y, z = r3c1
    if x + y > 10:
        next_x = (x * 3 - y) % 12
        next_y = (y * 2 + 4) % 12
        next_z = (z + x) % 12
    else: # x + y <= 10
        next_x = (x * 2 + y) % 12
        next_y = (y * 3 - 2) % 12
        next_z = (z * 2) % 12
    
    r3c2 = [next_x, next_y, next_z]
    print(f"Calculated middle triplet of Row 3: {r3c2}")

    # --- Calculate the rightmost triplet of Row 3 (M[2][2]) ---
    x, y, z = r3c2
    if x + y > 10:
        next_x = (x * 3 - y) % 12
        next_y = (y * 2 + 4) % 12
        next_z = (z + x) % 12
    else: # x + y <= 10
        next_x = (x * 2 + y) % 12
        next_y = (y * 3 - 2) % 12
        next_z = (z * 2) % 12

    r3c3 = [next_x, next_y, next_z]
    print(f"Calculated rightmost triplet of Row 3: {r3c3}")

    # --- Sum the elements of the missing triplets ---
    total_sum = sum(r3c2) + sum(r3c3)

    # --- Display the final sum and equation ---
    print("\nCalculating the sum of the missing elements:")
    sum_equation = f"{r3c2[0]} + {r3c2[1]} + {r3c2[2]} + {r3c3[0]} + {r3c3[1]} + {r3c3[2]} = {total_sum}"
    print(sum_equation)

solve_puzzle()
<<<24>>>