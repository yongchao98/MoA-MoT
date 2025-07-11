import math

def solve_cutting_problem():
    """
    Calculates the maximum value based on the user's problem formulation.

    The formulation has flaws, particularly that T1 and B2 pieces cannot coexist,
    and T1 pieces are packed very inefficiently. The optimal solution under the
    given rules involves maximizing B2 pieces and filling the rest with B1s.
    """

    # Prices for each product
    price_b2 = 150
    price_b1 = 1
    price_t1 = 5

    # --- Step 1: Calculate the maximum number of B2 balls ---
    # Billet dimensions in 0.5cm grid units: 32x22x8
    # B2 (radius 4) must be centered at z=4.
    # B2-B2 constraint: center-to-center distance >= 8.
    # We can place B2s on a grid with spacing 8.
    
    # x-dimension: Billet is 32 units wide. B2 center range is [4, 28].
    # Possible x centers: 4, 12, 20, 28. (4 positions)
    num_b2_x = math.floor((28 - 4) / 8) + 1

    # y-dimension: Billet is 22 units deep. B2 center range is [4, 18].
    # Possible y centers: 4, 12. (2 positions)
    num_b2_y = math.floor((18 - 4) / 8) + 1

    # z-dimension: Billet is 8 units high. B2 center range is z=4. (1 position)
    num_b2_z = 1
    
    num_b2 = num_b2_x * num_b2_y * num_b2_z
    
    # --- Step 2: Calculate B1 balls in remaining space ---
    # The 8 B2s occupy a volume of 16x8x4 cm (32x16x8 grid units).
    # This leaves a main rectangular block of 16x3x4 cm (32x6x8 grid units).
    # In grid coordinates, this block is x:[0,32], y:[16,22], z:[0,8].

    # B1 balls (radius 1) can be packed on a grid with spacing 2.
    # B1-B2 constraint is met for any B1 in this block, as the y-separation is sufficient.
    
    # B1s in the rectangular block:
    # x-centers in [1,31] with spacing 2: (31-1)/2 + 1 = 16
    b1_in_block_x = math.floor((31 - 1) / 2) + 1
    # y-centers in [17,21] with spacing 2: 17, 19, 21 => 3
    b1_in_block_y = math.floor((21 - 17) / 2) + 1
    # z-centers in [1,7] with spacing 2: 1, 3, 5, 7 => 4
    b1_in_block_z = math.floor((7 - 1) / 2) + 1
    
    num_b1_in_block = b1_in_block_x * b1_in_block_y * b1_in_block_z

    # B1s in the interstices between B2 balls:
    # The B2s are on a 4x2 grid in the x-y plane. This creates 3x1 = 3 interstitial columns.
    # (e.g., centered at x=8, y=8; x=16, y=8; x=24, y=8)
    num_interstitial_xy = (num_b2_x - 1) * (num_b2_y - 1)
    
    # At each interstitial location, we can stack B1s vertically.
    # The z-centers can be 1, 3, 5, 7 => 4 positions.
    num_interstitial_z = b1_in_block_z
    
    num_b1_in_interstices = num_interstitial_xy * num_interstitial_z
    
    # Total B1 count
    num_b1 = num_b1_in_block + num_b1_in_interstices

    # --- Step 3: T1 pieces ---
    # As reasoned in the plan, no T1 pieces can be included in the optimal solution
    # under the given constraints.
    num_t1 = 0
    
    # --- Step 4: Calculate final value and print the equation ---
    total_value = (num_b2 * price_b2) + (num_b1 * price_b1) + (num_t1 * price_t1)
    
    print("The problem formulation has flaws, making a mix of T1 and B2 pieces impossible.")
    print("Based on the provided rules, the highest value is achieved with a mix of B2 and B1 pieces.")
    print("\nCalculation for the highest value:")
    print(f"{num_b2} * {price_b2} (B2) + {num_b1} * {price_b1} (B1) + {num_t1} * {price_t1} (T1) = {total_value}")
    print(f"{num_b2 * price_b2} + {num_b1 * price_b1} + {num_t1 * price_t1} = {total_value}")
    
    # Final answer in the required format
    print("\nThe highest valid solution is:")
    print(f"<<<{total_value}>>>")

solve_cutting_problem()