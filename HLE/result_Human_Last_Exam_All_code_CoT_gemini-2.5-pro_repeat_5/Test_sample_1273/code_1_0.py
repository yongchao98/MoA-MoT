import math

def solve_billet_problem():
    """
    Calculates the maximum value that can be obtained by cutting a steel billet
    into specified products.
    """
    # Billet dimensions in cm
    billet_length = 8
    billet_width = 8
    billet_height = 4

    # Product 1: B2 Ball
    b2_radius = 2
    b2_price = 150
    # A ball with 2cm radius has a 4cm diameter, requiring a 4x4x4 cm block
    b2_block_dim = b2_radius * 2

    # The problem is to maximize value. Let's analyze the most valuable item, B2.
    # We will pursue a strategy of maximizing the number of B2 balls.

    # Calculate how many B2 balls can be cut from the billet
    num_b2_length = math.floor(billet_length / b2_block_dim)
    num_b2_width = math.floor(billet_width / b2_block_dim)
    num_b2_height = math.floor(billet_height / b2_block_dim)
    
    total_b2 = num_b2_length * num_b2_width * num_b2_height

    # This arrangement of B2 blocks (2 along length, 2 along width, 1 along height)
    # uses the entire 8x8x4 billet, leaving no space for other items.
    # Therefore, the maximum value is obtained by producing only B2 balls.

    max_value = total_b2 * b2_price

    print("Optimal Strategy: Maximize the production of the most valuable product, B2 balls.")
    print(f"Billet Dimensions: {billet_length}x{billet_width}x{billet_height} cm")
    print(f"B2 Ball Requirements: Needs a {b2_block_dim}x{b2_block_dim}x{b2_block_dim} cm block of steel.")
    print("\nCalculation of an optimal solution:")
    print(f"Number of B2 balls along length = floor({billet_length} / {b2_block_dim}) = {num_b2_length}")
    print(f"Number of B2 balls along width = floor({billet_width} / {b2_block_dim}) = {num_b2_width}")
    print(f"Number of B2 balls along height = floor({billet_height} / {b2_block_dim}) = {num_b2_height}")
    print(f"Total number of B2 balls = {num_b2_length} * {num_b2_width} * {num_b2_height} = {total_b2}")
    
    print("\nFinal Equation for Maximum Value:")
    print(f"Max Value = (Total B2 balls) * (Price per B2 ball)")
    print(f"Max Value = {total_b2} * {b2_price} = {max_value}")

solve_billet_problem()