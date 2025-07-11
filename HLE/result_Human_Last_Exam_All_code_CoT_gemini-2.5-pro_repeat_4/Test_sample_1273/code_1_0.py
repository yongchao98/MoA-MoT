def solve():
    """
    Calculates the maximum value that can be obtained by cutting a steel billet
    into specified products.
    """

    # Billet dimensions
    billet_dims = (8, 8, 4)

    # Product specifications: (name, required_cube_side, price)
    products = {
        'B2': {'side': 4.0, 'price': 150},
        'B1': {'side': 1.0, 'price': 1},
        'T1': {'side': 0.8, 'price': 1}
    }

    # --- Step 1: Calculate Value Densities to find the best product ---
    # Value Density = Price / Volume of the required cube
    # B2 Density: 150 / (4*4*4) = 2.34375
    # B1 Density: 1 / (1*1*1) = 1.0
    # T1 Density: 1 / (0.8*0.8*0.8) = 1.953125
    # B2 is the most value-dense product. We should prioritize it.

    # --- Step 2: Calculate how many B2 balls can be made ---
    b2_side = products['B2']['side']
    num_b2_x = int(billet_dims[0] / b2_side)
    num_b2_y = int(billet_dims[1] / b2_side)
    num_b2_z = int(billet_dims[2] / b2_side)
    
    total_b2 = num_b2_x * num_b2_y * num_b2_z
    
    # --- Step 3: Check for remaining material ---
    # Volume of billet = 8 * 8 * 4 = 256
    # Volume used for B2 = 4 * (4 * 4 * 4) = 4 * 64 = 256
    # No material is left.

    # --- Step 4: Calculate the maximum total value ---
    b2_price = products['B2']['price']
    max_value = total_b2 * b2_price

    print("The optimal strategy is to produce as many B2 balls as possible.")
    print(f"The billet (8x8x4 cm) can be cut into {total_b2} cubes of size 4x4x4 cm.")
    print("Each cube is used to make one B2 ball, with no material left over.")
    print("The calculation for the maximum value is:")
    
    # --- Final Output ---
    # The prompt requires printing each number in the final equation.
    print(f"{total_b2} * {b2_price} = {max_value}")

solve()