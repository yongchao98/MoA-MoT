import math

def solve_billet_problem():
    """
    Calculates the maximum value obtainable from cutting a steel billet,
    assuming leftover material is fully reusable.
    """
    # Billet dimensions
    billet_l, billet_w, billet_h = 8, 8, 4
    billet_volume = billet_l * billet_w * billet_h
    
    # Product specifications
    # B2: Ball, radius 2cm, price 150
    b2_radius = 2
    b2_diameter = b2_radius * 2
    b2_price = 150
    b2_volume = (4/3) * math.pi * (b2_radius**3)

    # B1: Ball, diameter 1cm, price 1
    b1_radius = 0.5
    b1_price = 1
    b1_volume = (4/3) * math.pi * (b1_radius**3)
    
    # T1: Cube, side 0.8cm, price 1
    t1_side = 0.8
    t1_price = 1
    t1_volume = t1_side**3

    # Step 1: Maximize the most valuable product, B2, based on physical fit.
    # A B2 ball requires a 4x4x4 cm space.
    num_b2 = math.floor(billet_l / b2_diameter) * math.floor(billet_w / b2_diameter) * math.floor(billet_h / b2_diameter)
    value_from_b2 = num_b2 * b2_price
    volume_used_by_b2 = num_b2 * b2_volume

    # Step 2: Use the remaining volume for the next most value-efficient product.
    # Value densities (value per cm^3 of product material):
    # B2: 150 / 33.51 ~= 4.47
    # T1: 1 / 0.512 ~= 1.95
    # B1: 1 / 0.523 ~= 1.91
    # T1 is more efficient than B1.
    remaining_volume = billet_volume - volume_used_by_b2
    num_t1 = math.floor(remaining_volume / t1_volume)
    value_from_t1 = num_t1 * t1_price
    
    # Step 3: Calculate total value
    total_value = value_from_b2 + value_from_t1

    print(f"Strategy: Maximize B2 balls, then fill remaining volume with T1 cubes.")
    print(f"1. Number of B2 balls that can be cut: {num_b2}")
    print(f"   Value from B2 balls: {num_b2} * {b2_price} = {value_from_b2}")
    print(f"2. Volume of billet: {billet_volume} cm^3")
    print(f"   Volume of {num_b2} B2 balls: {volume_used_by_b2:.2f} cm^3")
    print(f"   Remaining volume: {remaining_volume:.2f} cm^3")
    print(f"3. Number of T1 cubes from remaining volume: {num_t1}")
    print(f"   Value from T1 cubes: {num_t1} * {t1_price} = {value_from_t1}")
    print("\nFinal Calculation:")
    print(f"Total Value = (Value from B2) + (Value from T1)")
    print(f"Total Value = {value_from_b2} + {value_from_t1} = {total_value}")
    
    # Final equation as requested
    print("\nFinal equation with integer numbers:")
    print(f"{num_b2} * {b2_price} + {num_t1} * {t1_price} = {int(total_value)}")

solve_billet_problem()