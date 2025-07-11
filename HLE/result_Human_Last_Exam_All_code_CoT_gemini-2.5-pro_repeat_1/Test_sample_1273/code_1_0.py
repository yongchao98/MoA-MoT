import math

def solve():
    """
    Calculates the maximum value that can be obtained by cutting a steel billet
    into specified products.
    """
    # Billet dimensions
    billet_dims = (8, 8, 4)

    # Product specifications: (name, required_block_dims, price)
    products = {
        'B2': {'dims': (4, 4, 4), 'price': 150},
        'B1': {'dims': (1, 1, 1), 'price': 1},
        'T1': {'dims': (0.8, 0.8, 0.8), 'price': 1}
    }

    # Helper function to calculate max value from a given block using small items
    def get_max_value_from_block(block_dims):
        if not all(d > 0 for d in block_dims):
            return 0, 0, 0 # value, num_b1, num_t1

        # Calculate how many B1 balls can be cut
        b1_dims = products['B1']['dims']
        num_b1 = (math.floor(block_dims[0] / b1_dims[0]) *
                  math.floor(block_dims[1] / b1_dims[1]) *
                  math.floor(block_dims[2] / b1_dims[2]))
        value_b1 = num_b1 * products['B1']['price']

        # Calculate how many T1 cubes can be cut
        t1_dims = products['T1']['dims']
        num_t1 = (math.floor(block_dims[0] / t1_dims[0]) *
                  math.floor(block_dims[1] / t1_dims[1]) *
                  math.floor(block_dims[2] / t1_dims[2]))
        value_t1 = num_t1 * products['T1']['price']
        
        # Return the max value and the corresponding item counts
        if value_t1 > value_b1:
            return value_t1, 0, num_t1
        else:
            return value_b1, num_b1, 0

    # Main logic
    b2_dims = products['B2']['dims']
    max_b2 = (math.floor(billet_dims[0] / b2_dims[0]) *
              math.floor(billet_dims[1] / b2_dims[1]) *
              math.floor(billet_dims[2] / b2_dims[2]))

    best_strategy = {
        'total_value': 0,
        'num_b2': 0,
        'num_b1': 0,
        'num_t1': 0
    }

    print("Analyzing cutting strategies:")

    # Iterate through all possible numbers of B2 balls
    for num_b2 in range(max_b2, -1, -1):
        value_b2 = num_b2 * products['B2']['price']

        # The remaining volume can be considered as (max_b2 - num_b2) blocks of 4x4x4
        # This is a simplification, but valid since the billet perfectly tiles into 4x4x4 blocks.
        num_remaining_blocks = max_b2 - num_b2
        
        rem_value = 0
        rem_num_b1 = 0
        rem_num_t1 = 0

        if num_remaining_blocks > 0:
            # The remainder can be cut from one large block or multiple smaller ones.
            # Since T1s pack perfectly into a 4x4x4 block (5x5x5=125), we can calculate
            # the value from a single 4x4x4 block and multiply.
            single_block_value, _, single_block_num_t1 = get_max_value_from_block(b2_dims)
            rem_value = num_remaining_blocks * single_block_value
            rem_num_t1 = num_remaining_blocks * single_block_num_t1

        total_value = value_b2 + rem_value
        
        print(f"- Trying {num_b2} B2 ball(s):")
        print(f"  - Value from B2 balls: {num_b2} * {products['B2']['price']} = {value_b2}")
        print(f"  - Leftover material is equivalent to {num_remaining_blocks} block(s) of size {b2_dims}.")
        print(f"  - Max value from leftover material (using T1 cubes): {rem_value}")
        print(f"  - Total potential value: {value_b2} + {rem_value} = {total_value}\n")

        if total_value > best_strategy['total_value']:
            best_strategy['total_value'] = total_value
            best_strategy['num_b2'] = num_b2
            best_strategy['num_b1'] = rem_num_b1
            best_strategy['num_t1'] = rem_num_t1

    print("---")
    print("Best strategy found:")
    print(f"B2 Balls (2cm radius): {best_strategy['num_b2']}")
    print(f"B1 Balls (1cm diameter): {best_strategy['num_b1']}")
    print(f"T1 Cubes (0.8cm side): {best_strategy['num_t1']}")
    print("\nFinal Equation:")
    print(f"{best_strategy['num_b2']} * {products['B2']['price']} + "
          f"{best_strategy['num_b1']} * {products['B1']['price']} + "
          f"{best_strategy['num_t1']} * {products['T1']['price']} = "
          f"{best_strategy['total_value']}")

solve()