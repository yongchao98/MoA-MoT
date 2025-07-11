import math

def solve():
    """
    Calculates the maximum value that can be obtained by cutting a steel billet
    according to the problem's specifications.
    """
    # Billet dimensions
    billet_dims = (8, 8, 4)
    billet_volume = billet_dims[0] * billet_dims[1] * billet_dims[2]

    # Product specifications
    # B2: Ball, 2cm radius, price 150
    b2_radius = 2
    b2_diameter = b2_radius * 2
    b2_box_dims = (b2_diameter, b2_diameter, b2_diameter)
    b2_price = 150
    b2_volume = (4/3) * math.pi * (b2_radius**3)

    # T1: Cube, 0.8cm side, price 1
    t1_side = 0.8
    t1_price = 1
    t1_volume = t1_side**3

    # --- Strategy: Maximize B2 balls first, then use waste for T1 cubes ---

    # 1. Calculate how many B2 balls can be cut from the billet
    num_b2_x = math.floor(billet_dims[0] / b2_box_dims[0])
    num_b2_y = math.floor(billet_dims[1] / b2_box_dims[1])
    num_b2_z = math.floor(billet_dims[2] / b2_box_dims[2])
    num_b2 = num_b2_x * num_b2_y * num_b2_z

    # 2. Calculate the value from these B2 balls
    value_from_b2 = num_b2 * b2_price

    # 3. Calculate the waste material volume
    # The entire billet is used for the bounding boxes of the 4 B2 balls.
    # Waste = Volume of Billet - Volume of 4 spheres
    total_b2_volume = num_b2 * b2_volume
    waste_volume = billet_volume - total_b2_volume

    # 4. Calculate how many T1 cubes can be made from the waste
    num_t1_from_waste = math.floor(waste_volume / t1_volume)

    # 5. Calculate the value from these T1 cubes
    value_from_t1 = num_t1_from_waste * t1_price

    # 6. Calculate the total maximum value
    total_value = value_from_b2 + value_from_t1

    # Print the breakdown of the calculation for the final answer
    print("--- Optimal Strategy Calculation ---")
    print(f"Number of B2 balls (2cm radius): {num_b2}")
    print(f"Number of T1 cubes (0.8cm side) from waste: {num_t1_from_waste}")
    print("\n--- Final Equation ---")
    print(f"{num_b2} * {b2_price} (B2) + {num_t1_from_waste} * {t1_price} (T1) = {value_from_b2} + {value_from_t1} = {total_value}")

solve()