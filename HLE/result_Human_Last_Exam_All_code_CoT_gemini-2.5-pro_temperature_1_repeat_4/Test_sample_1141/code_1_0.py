import math

def solve():
    """
    Solves the billet cutting optimization problem based on the user's provided constraints.
    """

    # --- Product Data ---
    prices = {'B2': 150, 'B1': 1, 'T1': 5}
    
    # --- Analysis for each pure strategy ---

    # 1. B2 (Balls, 2cm radius) Only
    # Billet is 16x11x4 cm. B2 diameter is 4cm.
    # Center must be at z=2cm, so it's a 2D packing problem in a 16x11 rectangle.
    # We can fit 4 balls in a row along the 16cm side.
    # A staggered packing allows for 3 rows (e.g., at y=2, y=5.5, y=9)
    # with 4, 3, and 4 balls respectively.
    num_b2 = 11
    value_b2 = num_b2 * prices['B2']
    
    # 2. B1 (Balls, 1cm diameter) Only
    # B1 diameter is 1cm. We can pack them in a grid.
    num_b1_x = math.floor(16 / 1)
    num_b1_y = math.floor(11 / 1)
    num_b1_z = math.floor(4 / 1)
    num_b1 = num_b1_x * num_b1_y * num_b1_z
    value_b1 = num_b1 * prices['B1']

    # 3. T1 (Cubes, 1cm side) Only
    # User constraint: min(|xi-xj|, |yi-yj|, |zi-zj|) >= 2 (in 0.5cm grid units)
    # This means centers must be separated by >= 2 units on ALL axes.
    # We find the max number of points we can place in the space following this rule.
    # The number of such points is min(num_x_positions, num_y_positions, num_z_positions)
    # Num x positions: floor((31-1)/2) + 1 = 16
    # Num y positions: floor((21-1)/2) + 1 = 11
    # Num z positions: floor((7-1)/2) + 1 = 4
    num_t1 = min(16, 11, 4)
    value_t1 = num_t1 * prices['T1']

    # --- Mixed Strategy Analysis ---
    # As analyzed in the text, mixing strategies is not optimal under the given rules.
    # B2s take the full height, preventing any other parts.
    
    # --- Conclusion ---
    max_value = max(value_b2, value_b1, value_t1)

    print("The problem formulation for T1 cubes is geometrically incorrect for standard cutting.")
    print("However, solving the problem strictly by the rules provided:")
    print(f"- Pure B2 strategy yields: {value_b2}")
    print(f"- Pure B1 strategy yields: {value_b1}")
    print(f"- Pure T1 strategy yields: {value_t1}")
    print("\nThe optimal strategy is to produce only B2 balls.")
    
    # Generate the equation string
    sum_parts = [str(prices['B2'])] * num_b2
    equation = " + ".join(sum_parts)
    
    print("\nThe final equation for the highest value is:")
    print(f"{equation} = {max_value}")

solve()
<<<1650>>>