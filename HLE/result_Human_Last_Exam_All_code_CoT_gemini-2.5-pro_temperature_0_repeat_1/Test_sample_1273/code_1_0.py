import math

def solve():
    """
    Calculates the maximum value obtainable by cutting a steel billet
    into specified products.
    """

    # Billet dimensions
    billet_dims = (8, 8, 4)

    # Product specifications: (name, price, required_cube_dims)
    products = {
        'B2': {'price': 150, 'dims': (4, 4, 4)}, # Ball, 2cm radius -> 4cm diameter
        'B1': {'price': 1, 'dims': (1, 1, 1)},   # Ball, 1cm diameter
        'T1': {'price': 1, 'dims': (0.8, 0.8, 0.8)} # Cube, 0.8cm side
    }

    # --- Analysis ---
    # Comparing B1 and T1:
    # For a given volume, T1 is more profitable than B1 because it's smaller
    # for the same price and all billet dimensions are perfect multiples of
    # both 1cm and 0.8cm, meaning no packing advantage for B1.
    # Value per cm^3 of required space:
    # B1: 1 / (1*1*1) = 1
    # T1: 1 / (0.8*0.8*0.8) = 1 / 0.512 = ~1.95
    # So, we prioritize T1 over B1 for any leftover space.

    # --- Strategies ---
    # We will test strategies based on the number of B2 balls, as they are the
    # most valuable item. We fill the rest of the space with T1 cubes.

    strategies = {}

    # Strategy 1: Maximize B2 balls
    num_b2_s1 = (math.floor(billet_dims[0] / products['B2']['dims'][0]) *
                 math.floor(billet_dims[1] / products['B2']['dims'][1]) *
                 math.floor(billet_dims[2] / products['B2']['dims'][2]))
    value_s1 = num_b2_s1 * products['B2']['price']
    strategies['4 B2 balls'] = {'value': value_s1, 'num_b2': num_b2_s1, 'num_t1': 0}

    # Strategy 2: 3 B2 balls, rest with T1
    num_b2_s2 = 3
    value_b2_s2 = num_b2_s2 * products['B2']['price']
    # Remaining volume is a 4x4x4 block
    rem_dims_s2 = (4, 4, 4)
    num_t1_s2 = (math.floor(rem_dims_s2[0] / products['T1']['dims'][0]) *
                 math.floor(rem_dims_s2[1] / products['T1']['dims'][1]) *
                 math.floor(rem_dims_s2[2] / products['T1']['dims'][2]))
    value_t1_s2 = num_t1_s2 * products['T1']['price']
    strategies['3 B2 balls + T1 cubes'] = {'value': value_b2_s2 + value_t1_s2, 'num_b2': num_b2_s2, 'num_t1': num_t1_s2}

    # Strategy 3: 2 B2 balls, rest with T1
    num_b2_s3 = 2
    value_b2_s3 = num_b2_s3 * products['B2']['price']
    # Remaining volume is an 8x4x4 block
    rem_dims_s3 = (8, 4, 4)
    num_t1_s3 = (math.floor(rem_dims_s3[0] / products['T1']['dims'][0]) *
                 math.floor(rem_dims_s3[1] / products['T1']['dims'][1]) *
                 math.floor(rem_dims_s3[2] / products['T1']['dims'][2]))
    value_t1_s3 = num_t1_s3 * products['T1']['price']
    strategies['2 B2 balls + T1 cubes'] = {'value': value_b2_s3 + value_t1_s3, 'num_b2': num_b2_s3, 'num_t1': num_t1_s3}

    # Strategy 4: 1 B2 ball, rest with T1
    num_b2_s4 = 1
    value_b2_s4 = num_b2_s4 * products['B2']['price']
    # Remaining volume is an L-shape, which can be cut from two blocks: 4x8x4 and 4x4x4
    rem_dims1_s4 = (4, 8, 4)
    rem_dims2_s4 = (4, 4, 4)
    num_t1_s4_1 = (math.floor(rem_dims1_s4[0] / products['T1']['dims'][0]) *
                   math.floor(rem_dims1_s4[1] / products['T1']['dims'][1]) *
                   math.floor(rem_dims1_s4[2] / products['T1']['dims'][2]))
    num_t1_s4_2 = (math.floor(rem_dims2_s4[0] / products['T1']['dims'][0]) *
                   math.floor(rem_dims2_s4[1] / products['T1']['dims'][1]) *
                   math.floor(rem_dims2_s4[2] / products['T1']['dims'][2]))
    num_t1_s4 = num_t1_s4_1 + num_t1_s4_2
    value_t1_s4 = num_t1_s4 * products['T1']['price']
    strategies['1 B2 ball + T1 cubes'] = {'value': value_b2_s4 + value_t1_s4, 'num_b2': num_b2_s4, 'num_t1': num_t1_s4}

    # Strategy 5: 0 B2 balls, all T1
    num_b2_s5 = 0
    num_t1_s5 = (math.floor(billet_dims[0] / products['T1']['dims'][0]) *
                 math.floor(billet_dims[1] / products['T1']['dims'][1]) *
                 math.floor(billet_dims[2] / products['T1']['dims'][2]))
    value_t1_s5 = num_t1_s5 * products['T1']['price']
    strategies['All T1 cubes'] = {'value': value_t1_s5, 'num_b2': num_b2_s5, 'num_t1': num_t1_s5}

    # --- Find the best strategy ---
    best_strategy_name = max(strategies, key=lambda k: strategies[k]['value'])
    best_strategy_details = strategies[best_strategy_name]
    max_value = best_strategy_details['value']
    
    print("The best strategy is to cut the billet to produce the following:")
    
    final_equation = []
    if best_strategy_details['num_b2'] > 0:
        b2_value = best_strategy_details['num_b2'] * products['B2']['price']
        print(f"- {best_strategy_details['num_b2']} B2 balls")
        final_equation.append(f"{best_strategy_details['num_b2']} * {products['B2']['price']}")

    if best_strategy_details.get('num_t1', 0) > 0:
        t1_value = best_strategy_details['num_t1'] * products['T1']['price']
        print(f"- {best_strategy_details['num_t1']} T1 cubes")
        final_equation.append(f"{best_strategy_details['num_t1']} * {products['T1']['price']}")

    print("\nThe calculation for the total value is:")
    print(f"{' + '.join(final_equation)} = {int(max_value)}")

solve()