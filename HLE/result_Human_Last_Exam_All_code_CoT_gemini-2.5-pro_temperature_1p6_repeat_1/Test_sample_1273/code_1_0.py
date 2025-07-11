import math

def solve_cutting_problem():
    """
    Analyzes different cutting strategies for a steel billet to maximize product value.
    """
    # --- Define problem parameters ---
    billet_dims = {'length': 8, 'width': 8, 'height': 4}
    
    products = {
        'B2': {'name': 'Ball (2cm radius)', 'stock_dims': (4, 4, 4), 'price': 150},
        'B1': {'name': 'Ball (1cm diameter)', 'stock_dims': (1, 1, 1), 'price': 1},
        'T1': {'name': 'Cube (0.8cm side)', 'stock_dims': (0.8, 0.8, 0.8), 'price': 1}
    }

    print("Analyzing cutting options for an 8x8x4 cm billet.")
    print("-" * 50)

    # --- Strategy 1: Cut only the most valuable item, B2 balls ---
    b2 = products['B2']
    num_b2_x = math.floor(billet_dims['length'] / b2['stock_dims'][0])
    num_b2_y = math.floor(billet_dims['width'] / b2['stock_dims'][1])
    num_b2_z = math.floor(billet_dims['height'] / b2['stock_dims'][2])
    total_b2 = num_b2_x * num_b2_y * num_b2_z
    value_b2_only = total_b2 * b2['price']
    
    print("Strategy 1: Fill the entire billet with B2 balls.")
    print(f"A B2 ball requires a {b2['stock_dims'][0]}x{b2['stock_dims'][1]}x{b2['stock_dims'][2]} cm block.")
    print(f"Number of B2 balls = floor(8/4) * floor(8/4) * floor(4/4) = {num_b2_x} * {num_b2_y} * {num_b2_z} = {total_b2}")
    print(f"This strategy uses the entire billet.")
    print(f"Total Value = {total_b2} * {b2['price']} = {value_b2_only}")
    print("-" * 50)
    
    max_value = value_b2_only
    best_strategy_desc = f"Cut {total_b2} B2 balls."
    final_equation_nums = [total_b2, b2['price'], value_b2_only]

    # --- Strategy 2: Cut only T1 cubes ---
    t1 = products['T1']
    num_t1_x = math.floor(billet_dims['length'] / t1['stock_dims'][0])
    num_t1_y = math.floor(billet_dims['width'] / t1['stock_dims'][1])
    num_t1_z = math.floor(billet_dims['height'] / t1['stock_dims'][2])
    total_t1 = num_t1_x * num_t1_y * num_t1_z
    value_t1_only = total_t1 * t1['price']

    print("Strategy 2: Fill the entire billet with T1 cubes.")
    print(f"Number of T1 cubes = floor(8/0.8) * floor(8/0.8) * floor(4/0.8) = {num_t1_x} * {num_t1_y} * {num_t1_z} = {total_t1}")
    print(f"Total Value = {total_t1} * {t1['price']} = {value_t1_only}")
    print("-" * 50)
    
    # --- Strategy 3: Mixed - 3 B2s + T1 cubes ---
    print("Strategy 3: Cut 3 B2 balls and fill the remainder with T1 cubes.")
    value_from_3_b2 = 3 * b2['price']
    # The remainder is one 4x4x4 block
    rem_block_dims = (4, 4, 4)
    rem_t1_x = math.floor(rem_block_dims[0] / t1['stock_dims'][0])
    rem_t1_y = math.floor(rem_block_dims[1] / t1['stock_dims'][1])
    rem_t1_z = math.floor(rem_block_dims[2] / t1['stock_dims'][2])
    rem_total_t1 = rem_t1_x * rem_t1_y * rem_t1_z
    rem_value_t1 = rem_total_t1 * t1['price']
    value_mixed_3_b2 = value_from_3_b2 + rem_value_t1
    print(f"Value from 3 B2 balls = 3 * 150 = {value_from_3_b2}")
    print(f"The remainder is a 4x4x4 cm block, which can fit {rem_total_t1} T1 cubes (Value = {rem_value_t1}).")
    print(f"Total Value = {value_from_3_b2} + {rem_value_t1} = {value_mixed_3_b2}")
    print("-" * 50)

    # --- Final Conclusion ---
    print("Conclusion:")
    print("Comparing the calculated values (600, 500, 575, etc.), the maximum achievable value is 600.")
    print(f"The best strategy is to {best_strategy_desc}")
    print("The final calculation for the best option is:")
    # Using the requirement to print each number in the final equation
    print(f"{final_equation_nums[0]} * {final_equation_nums[1]} = {final_equation_nums[2]}")

solve_cutting_problem()
<<<A>>>