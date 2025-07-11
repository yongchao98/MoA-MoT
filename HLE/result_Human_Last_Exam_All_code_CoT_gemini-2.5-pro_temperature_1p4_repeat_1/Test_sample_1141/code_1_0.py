import math

def solve_cutting_problem():
    """
    Solves the billet cutting problem by finding the product with the highest
    value density and calculating the maximum value achievable.
    """
    # 1. Define billet and product properties
    billet_dims = {'length': 16, 'width': 11, 'height': 4}
    products = {
        'B2': {'type': 'ball', 'radius': 2, 'price': 150},
        'B1': {'type': 'ball', 'radius': 0.5, 'price': 1}, # 1cm diameter
        'T1': {'type': 'cube', 'side': 1, 'price': 5}
    }

    # 2. Calculate volume and value density for each product
    for name, props in products.items():
        if props['type'] == 'ball':
            volume = (4/3) * math.pi * (props['radius'] ** 3)
        elif props['type'] == 'cube':
            volume = props['side'] ** 3
        
        props['volume'] = volume
        props['value_density'] = props['price'] / volume

    # 3. Identify the product with the highest value density
    best_product_name = max(products, key=lambda p: products[p]['value_density'])
    best_product = products[best_product_name]

    # 4. Determine the optimal strategy. Since T1 cubes are most value-dense
    # and can tile the space perfectly, filling the billet with them is optimal.
    if best_product['type'] == 'cube':
        # 5. Calculate total pieces and maximum value
        num_pieces_length = billet_dims['length'] // best_product['side']
        num_pieces_width = billet_dims['width'] // best_product['side']
        num_pieces_height = billet_dims['height'] // best_product['side']

        total_pieces = num_pieces_length * num_pieces_width * num_pieces_height
        price_per_piece = best_product['price']
        max_value = total_pieces * price_per_piece

        print("The optimal strategy is to cut the entire billet into T1 cubes.")
        print("The calculation for the maximum value is:")
        # Output each number in the final equation as requested
        print(f"{int(total_pieces)} * {int(price_per_piece)} = {int(max_value)}")
    else:
        # This case is not reached here, but is included for completeness
        print("The optimal product is not a cube, complex packing analysis would be required.")


solve_cutting_problem()