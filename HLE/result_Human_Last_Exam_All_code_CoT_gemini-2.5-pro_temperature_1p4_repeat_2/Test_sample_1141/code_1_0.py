import math

def solve_billet_packing():
    """
    Solves the user's specific cutting stock problem to maximize value.
    The function first places the most valuable items (B2 balls) and then
    greedily fills the remaining space with B1 balls.
    """
    
    # --- Item Properties ---
    # Prices are as specified by the user.
    B2_PRICE = 150
    B1_PRICE = 1
    T1_PRICE = 5

    # --- B2 Ball Placement ---
    # An efficient packing allows for 11 B2 balls.
    # All B2 balls are at z=4 as per the problem description.
    b2_centers = [
        (4, 4, 4), (12, 4, 4), (20, 4, 4), (28, 4, 4),
        (8, 11, 4), (16, 11, 4), (24, 11, 4),
        (4, 18, 4), (12, 18, 4), (20, 18, 4), (28, 18, 4)
    ]
    num_b2 = len(b2_centers)

    # --- Greedy B1 Ball Placement ---
    # Now, fill the remaining space with B1 balls.
    b1_centers = []
    
    # Pre-calculate squared distances for constraint checking to improve performance.
    B1_TO_B2_DIST_SQ = 25  # (radius_B1 + radius_B2)^2 = (1+4)^2 = 25
    B1_TO_B1_DIST_SQ = 4   # (radius_B1 + radius_B1)^2 = (1+1)^2 = 4
    
    # Define valid center ranges for B1 balls
    b1_x_range = (1, 31)
    b1_y_range = (1, 21)
    b1_z_range = (1, 7)
    
    # Iterate through all possible B1 center locations in a greedy fashion
    for z1 in range(b1_z_range[0], b1_z_range[1] + 1):
        if z1 == 4:  # Skip the z-level occupied by B2 balls
            continue
        for y1 in range(b1_y_range[0], b1_y_range[1] + 1):
            for x1 in range(b1_x_range[0], b1_x_range[1] + 1):
                
                is_valid_location = True
                
                # Check non-overlapping constraint with all B2 balls
                for xb2, yb2, zb2 in b2_centers:
                    dist_sq = (x1 - xb2)**2 + (y1 - yb2)**2 + (z1 - zb2)**2
                    if dist_sq < B1_TO_B2_DIST_SQ:
                        is_valid_location = False
                        break
                if not is_valid_location:
                    continue
                
                # Check non-overlapping constraint with previously placed B1 balls
                for xb1, yb1, zb1 in b1_centers:
                    dist_sq = (x1 - xb1)**2 + (y1 - yb1)**2 + (z1 - zb1)**2
                    if dist_sq < B1_TO_B1_DIST_SQ:
                        is_valid_location = False
                        break
                
                if is_valid_location:
                    b1_centers.append((x1, y1, z1))

    num_b1 = len(b1_centers)
    num_t1 = 0 # As established, no T1s can be placed with B2s

    # --- Calculate and Print Final Value ---
    total_value = (num_b2 * B2_PRICE) + (num_b1 * B1_PRICE) + (num_t1 * T1_PRICE)
    
    print("Based on the problem formulation, the highest value is achieved by combining B2 balls and B1 balls.")
    print(f"Number of B2 (2cm radius ball) pieces: {num_b2}")
    print(f"Number of B1 (1cm diameter ball) pieces: {num_b1}")
    print(f"Number of T1 (1cm side cube) pieces: {num_t1}")
    print("\nThe final calculation for the maximum value is:")
    print(f"{num_b2} * {B2_PRICE} + {num_b1} * {B1_PRICE} + {num_t1} * {T1_PRICE} = {total_value}")

if __name__ == '__main__':
    solve_billet_packing()