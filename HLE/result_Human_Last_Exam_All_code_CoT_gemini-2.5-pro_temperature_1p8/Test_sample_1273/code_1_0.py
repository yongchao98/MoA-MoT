import math

def solve():
    """
    Calculates the maximum value that can be obtained by cutting a steel billet
    into specified product shapes.
    """
    # Billet dimensions
    billet_dims = (8, 8, 4)

    # Product specifications: (name, price, required_block_dims)
    products = {
        'B2': {'price': 150, 'dims': (4, 4, 4)},
        'B1': {'price': 1, 'dims': (1, 1, 1)},
        'T1': {'price': 1, 'dims': (0.8, 0.8, 0.8)}
    }

    # As analyzed in the plan, T1 cubes are more space-efficient than B1 balls.
    # We will prioritize T1 for filling leftover space.
    filler_item = products['T1']

    def calculate_filler_value(block_dims):
        """Calculates how many filler items (T1) fit in a given block and their total value."""
        if not all(d > 0 for d in block_dims):
            return 0, 0
            
        num_x = math.floor(block_dims[0] / filler_item['dims'][0])
        num_y = math.floor(block_dims[1] / filler_item['dims'][1])
        num_z = math.floor(block_dims[2] / filler_item['dims'][2])
        
        num_items = num_x * num_y * num_z
        total_value = num_items * filler_item['price']
        return num_items, total_value

    max_b2_x = math.floor(billet_dims[0] / products['B2']['dims'][0])
    max_b2_y = math.floor(billet_dims[1] / products['B2']['dims'][1])
    max_b2_z = math.floor(billet_dims[2] / products['B2']['dims'][2])
    max_b2_possible = max_b2_x * max_b2_y * max_b2_z

    best_scenario = {
        'value': 0,
        'num_b2': 0,
        'num_t1': 0
    }
    
    print("Evaluating cutting scenarios:")
    # Loop through all possible numbers of B2 balls
    for num_b2 in range(max_b2_possible, -1, -1):
        value_b2 = num_b2 * products['B2']['price']
        
        # Calculate remaining space and the value that can be generated from it.
        # This requires analyzing how the blocks are cut.
        
        num_t1 = 0
        value_t1 = 0

        if num_b2 == 4:
            # 4 B2 blocks (2x2x1 layout) use the entire 8x8x4 billet.
            remaining_blocks = []
        elif num_b2 == 3:
            # 3 B2 blocks leave one 4x4x4 block.
            remaining_blocks = [(4, 4, 4)]
        elif num_b2 == 2:
            # 2 B2 blocks (e.g., side-by-side) leave an 8x4x4 block.
            remaining_blocks = [(8, 4, 4)]
        elif num_b2 == 1:
            # 1 B2 block leaves one 4x8x4 block and one 4x4x4 block.
            remaining_blocks = [(4, 8, 4), (4, 4, 4)]
        elif num_b2 == 0:
            # 0 B2 blocks leaves the entire billet.
            remaining_blocks = [billet_dims]
            
        for block in remaining_blocks:
            n, v = calculate_filler_value(block)
            num_t1 += n
            value_t1 += v

        total_value = value_b2 + value_t1
        
        print(f"- Scenario with {num_b2} B2s:")
        print(f"  - Value from B2s: {num_b2} * {products['B2']['price']} = {value_b2}")
        print(f"  - Value from T1s in leftover space: {num_t1} * {products['T1']['price']} = {value_t1}")
        print(f"  - Total Value: {total_value}\n")

        if total_value > best_scenario['value']:
            best_scenario['value'] = total_value
            best_scenario['num_b2'] = num_b2
            best_scenario['num_t1'] = num_t1

    print("-----------------------------------------")
    print("Final Result:")
    print(f"The maximum value achievable is {best_scenario['value']}.")
    print("This is achieved by the following cutting plan:")
    print(f"- Number of B2 balls (2cm radius): {best_scenario['num_b2']}")
    print(f"- Number of B1 balls (1cm diameter): 0")
    print(f"- Number of T1 cubes (0.8cm side): {best_scenario['num_t1']}")
    print("\nFinal equation for the best option:")
    print(f"{best_scenario['num_b2']} * {products['B2']['price']} + 0 * {products['B1']['price']} + {best_scenario['num_t1']} * {products['T1']['price']} = {best_scenario['value']}")

solve()
<<<A>>>