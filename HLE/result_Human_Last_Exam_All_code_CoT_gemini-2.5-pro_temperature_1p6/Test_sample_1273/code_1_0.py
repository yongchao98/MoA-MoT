def solve():
    """
    Calculates the maximum value that can be obtained by cutting a steel billet.
    """
    # Billet dimensions in cm
    billet_dim = (8, 8, 4)

    # Products: defined by their price and the dimensions of the cubic space required for cutting.
    products = {
        'B2': {'price': 150, 'dim': (4, 4, 4)},   # Ball with radius=2cm has diameter=4cm
        'B1': {'price': 1, 'dim': (1, 1, 1)},     # Ball with diameter=1cm
        'T1': {'price': 1, 'dim': (0.8, 0.8, 0.8)} # Cube with side=0.8cm
    }

    # Based on value density analysis, B2 is the most profitable product per unit of volume.
    # We will use a greedy strategy to fit as many B2 balls as possible.

    # Calculate how many B2 balls can be cut from the billet.
    b2_dim = products['B2']['dim']
    num_b2_x = int(billet_dim[0] / b2_dim[0])
    num_b2_y = int(billet_dim[1] / b2_dim[1])
    num_b2_z = int(billet_dim[2] / b2_dim[2])

    # Total number of B2 balls
    num_b2 = num_b2_x * num_b2_y * num_b2_z

    # The cubic blocks for the B2 balls fit perfectly into the billet,
    # leaving no remaining space for other products.
    # Therefore, the number of other products is 0.
    num_b1 = 0
    num_t1 = 0

    # Calculate the total maximum value
    total_value = (num_b2 * products['B2']['price']) + \
                  (num_b1 * products['B1']['price']) + \
                  (num_t1 * products['T1']['price'])

    # Print the final breakdown and the equation for the total value.
    print(f"The optimal strategy is to cut the billet into {num_b2} B2 balls.")
    print("This configuration utilizes the entire billet.")
    print("\nFinal Equation:")
    # The final printout includes each number in the equation as requested.
    print(f"{num_b2} * {products['B2']['price']} + {num_b1} * {products['B1']['price']} + {num_t1} * {products['T1']['price']} = {int(total_value)}")

solve()
<<<A>>>