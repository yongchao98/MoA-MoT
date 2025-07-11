import math

def solve_billet_optimization():
    """
    Calculates the maximum value that can be obtained from a steel billet
    by cutting it into three types of products.

    This solution assumes that all scrap material is reusable (volume conservation).
    """

    # Billet dimensions and volume
    billet_dims = (8, 8, 4)
    billet_volume = billet_dims[0] * billet_dims[1] * billet_dims[2]

    # Product specifications
    products = {
        'B2': {'type': 'ball', 'radius': 2, 'price': 150},
        'B1': {'type': 'ball', 'diameter': 1, 'price': 1},
        'T1': {'type': 'cube', 'side': 0.8, 'price': 1}
    }

    # Calculate volume and value density for each product
    for name, p in products.items():
        if p['type'] == 'ball':
            if 'radius' in p:
                radius = p['radius']
            else:
                radius = p['diameter'] / 2.0
            p['volume'] = (4/3) * math.pi * (radius ** 3)
        elif p['type'] == 'cube':
            p['volume'] = p['side'] ** 3
        p['value_density'] = p['price'] / p['volume']
        
    # --- Step 1: Maximize the most value-dense product (B2) ---
    # The number of B2 balls is limited by the billet's physical dimensions.
    # A B2 ball (radius 2) needs its center to be at least 2cm from any face.
    # Billet: 8x8x4. Center must be in [2, 6] x [2, 6] x [2, 2].
    # This defines a 4x4 plane for the centers.
    # The distance between centers of two r=2 spheres must be >= 4.
    # In a 4x4 plane, we can place 4 such centers.
    num_b2 = 4
    value_from_b2 = num_b2 * products['B2']['price']
    volume_used_by_b2 = num_b2 * products['B2']['volume']

    # --- Step 2: Use remaining volume for the next best product (T1) ---
    remaining_volume = billet_volume - volume_used_by_b2
    
    # The next most value-dense product is T1.
    num_t1 = math.floor(remaining_volume / products['T1']['volume'])
    value_from_t1 = num_t1 * products['T1']['price']
    
    # --- Step 3: Calculate total value ---
    total_value = value_from_b2 + value_from_t1

    # Output the final equation
    print("Maximization Plan:")
    print(f"- Prioritize B2 balls due to highest value density.")
    print(f"- The 8x8x4 cm billet can physically fit a maximum of {num_b2} B2 balls.")
    print(f"- These {num_b2} balls generate a value of {num_b2} * {products['B2']['price']} = {value_from_b2}.")
    print(f"- The steel used for B2 balls is {volume_used_by_b2:.2f} cm^3.")
    print(f"- The billet's total volume is {billet_volume} cm^3, leaving {remaining_volume:.2f} cm^3 of steel.")
    print(f"- This remaining steel is used for T1 cubes, the next most valuable product.")
    print(f"- Number of T1 cubes made: floor({remaining_volume:.2f} / {products['T1']['volume']}) = {num_t1}.")
    print(f"- These {num_t1} cubes generate a value of {num_t1} * {products['T1']['price']} = {value_from_t1}.")
    print("\nFinal Optimal Value Calculation:")
    print(f"{num_b2} * {products['B2']['price']} + {num_t1} * {products['T1']['price']} = {int(total_value)}")

solve_billet_optimization()