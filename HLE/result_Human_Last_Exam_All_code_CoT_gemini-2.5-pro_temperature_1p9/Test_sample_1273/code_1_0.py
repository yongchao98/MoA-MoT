import math

def solve_billet_problem():
    """
    Analyzes cutting strategies for a steel billet to maximize product value
    and prints the details of the best strategy.
    """

    # --- 1. Define Billet and Product Specifications ---
    billet_dims = {'length': 8, 'width': 8, 'height': 4}

    products = {
        'B2': {'name': 'B2 (Ball, 2cm radius)', 'block_dims': (4, 4, 4), 'price': 150},
        'B1': {'name': 'B1 (Ball, 1cm diameter)', 'block_dims': (1, 1, 1), 'price': 1},
        'T1': {'name': 'T1 (Cube, 0.8cm side)', 'block_dims': (0.8, 0.8, 0.8), 'price': 1}
    }

    # As determined in the plan, T1 is always better than B1 for filling space.
    # We will prioritize T1 for all leftover material.
    
    print("Analyzing cutting strategies to maximize value from an 8x8x4 cm billet...")
    print("-" * 30)

    best_scenario = {
        'value': 0,
        'num_b2': 0,
        'num_b1': 0,
        'num_t1': 0
    }

    # --- 2. Evaluate Strategies Based on Number of B2 Balls ---

    # Strategy 1: Cut maximum possible B2 balls (4)
    # 4 blocks of 4x4x4 fit perfectly into the 8x8x4 billet (2x2x1 layout).
    num_b2_s1 = math.floor(billet_dims['length'] / products['B2']['block_dims'][0]) * \
                math.floor(billet_dims['width'] / products['B2']['block_dims'][1]) * \
                math.floor(billet_dims['height'] / products['B2']['block_dims'][2])
    num_t1_s1 = 0
    value_s1 = num_b2_s1 * products['B2']['price']
    if value_s1 > best_scenario['value']:
        best_scenario.update({'value': value_s1, 'num_b2': num_b2_s1, 'num_t1': num_t1_s1})

    # Strategy 2: Cut 3 B2 balls
    # This leaves a single 4x4x4 cm block of remaining material.
    num_b2_s2 = 3
    rem_block_s2 = (4, 4, 4)
    num_t1_s2 = math.floor(rem_block_s2[0] / products['T1']['block_dims'][0]) * \
                math.floor(rem_block_s2[1] / products['T1']['block_dims'][1]) * \
                math.floor(rem_block_s2[2] / products['T1']['block_dims'][2])
    value_s2 = (num_b2_s2 * products['B2']['price']) + (num_t1_s2 * products['T1']['price'])
    if value_s2 > best_scenario['value']:
        best_scenario.update({'value': value_s2, 'num_b2': num_b2_s2, 'num_t1': num_t1_s2})
        
    # Strategy 3: Cut 2 B2 balls
    # Placing two 4x4x4 blocks side-by-side leaves a single 8x4x4 cm block.
    num_b2_s3 = 2
    rem_block_s3 = (8, 4, 4)
    num_t1_s3 = math.floor(rem_block_s3[0] / products['T1']['block_dims'][0]) * \
                math.floor(rem_block_s3[1] / products['T1']['block_dims'][1]) * \
                math.floor(rem_block_s3[2] / products['T1']['block_dims'][2])
    value_s3 = (num_b2_s3 * products['B2']['price']) + (num_t1_s3 * products['T1']['price'])
    if value_s3 > best_scenario['value']:
        best_scenario.update({'value': value_s3, 'num_b2': num_b2_s3, 'num_t1': num_t1_s3})

    # Strategy 4: Cut 0 B2 balls (fill entire billet with T1 cubes)
    num_b2_s4 = 0
    num_t1_s4 = math.floor(billet_dims['length'] / products['T1']['block_dims'][0]) * \
                math.floor(billet_dims['width'] / products['T1']['block_dims'][1]) * \
                math.floor(billet_dims['height'] / products['T1']['block_dims'][2])
    value_s4 = num_t1_s4 * products['T1']['price']
    if value_s4 > best_scenario['value']:
        best_scenario.update({'value': value_s4, 'num_b2': num_b2_s4, 'num_t1': num_t1_s4})

    # --- 3. Print the Best Found Scenario ---
    print("The best strategy found is:")
    print(f"  - Number of B2 balls: {best_scenario['num_b2']}")
    print(f"  - Number of B1 balls: {best_scenario['num_b1']} (less value-dense than T1)")
    print(f"  - Number of T1 cubes: {best_scenario['num_t1']}")
    
    print("\nThis strategy utilizes the entire billet to produce the most valuable items.")
    print("The final calculation for the total value is:")
    
    # Final equation printout
    b2_val = best_scenario['num_b2'] * products['B2']['price']
    b1_val = best_scenario['num_b1'] * products['B1']['price']
    t1_val = best_scenario['num_t1'] * products['T1']['price']
    
    print(
        f"({best_scenario['num_b2']} B2 * {products['B2']['price']}) + "
        f"({best_scenario['num_b1']} B1 * {products['B1']['price']}) + "
        f"({best_scenario['num_t1']} T1 * {products['T1']['price']}) = "
        f"{best_scenario['value']}"
    )

    print(f"{b2_val} + {b1_val} + {t1_val} = {best_scenario['value']}")
    
    print("\nComparing this to the options: A=600, B=640, C=648, D=838.")
    print("The maximum valid value is 600, which corresponds to option A.")


solve_billet_problem()