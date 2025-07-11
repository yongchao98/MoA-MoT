import numpy as np

def solve():
    """
    Calculates the maximal Thurston-Bennequin number for a given 5x5 grid diagram.
    """
    n = 5
    # The permutation sigma is defined by the positions of X's.
    # X's at (1,4), (2,5), (3,1), (4,2), (5,3)
    # sigma(i) = j means X in column i is in row j.
    # Using 0-based indexing for arrays, so subtract 1 from coordinates.
    # sigma = [sigma(1)-1, sigma(2)-1, ...]
    sigma = [3, 4, 0, 1, 2]

    # Calculate the inverse permutation sigma_inv
    sigma_inv = [0] * n
    for i in range(n):
        sigma_inv[sigma[i]] = i

    # --- Step 1: Calculate the writhe w(D) ---
    writhe = 0
    # A crossing at (i, j) exists if the vertical segment in column i
    # crosses the horizontal segment in row j.
    # Vertical segment i connects (i, i) to (i, sigma[i]).
    # Horizontal segment j connects (j, j) to (sigma_inv[j], j).
    for i in range(n):
        for j in range(n):
            # Check for crossing at (i, j) using 0-based indices
            min_y, max_y = min(i, sigma[i]), max(i, sigma[i])
            min_x, max_x = min(j, sigma_inv[j]), max(j, sigma_inv[j])
            
            # A crossing occurs if column i is between the horizontal segment's endpoints,
            # and row j is between the vertical segment's endpoints.
            if (min_y < j < max_y) and (min_x < i < max_x):
                # Vertical segment orientation
                v_sign = 1 if sigma[i] > i else -1
                # Horizontal segment orientation (from O at (j,j) to X at (sigma_inv[j],j))
                # So the direction is sgn(sigma_inv[j] - j)
                h_sign = 1 if sigma_inv[j] > j else -1
                
                # Sign of crossing is product of signs (vertical over horizontal)
                writhe += v_sign * h_sign

    # --- Step 2: Calculate the number of positive O-corners c_+(D) ---
    c_plus = 0
    for i in range(n):
        # O-corner at (i,i) is NE-type if sigma(i)>i and sigma_inv(i)>i
        is_NE = (sigma[i] > i) and (sigma_inv[i] > i)
        # O-corner at (i,i) is SW-type if sigma(i)<i and sigma_inv(i)<i
        is_SW = (sigma[i] < i) and (sigma_inv[i] < i)
        
        if is_NE or is_SW:
            c_plus += 1
            
    # --- Step 3: Calculate the Thurston-Bennequin number ---
    tb = writhe - c_plus
    
    print("Calculation Steps:")
    print(f"1. Writhe of the diagram, w(D) = {writhe}")
    print(f"2. Number of positive corners, c_+(D) = {c_plus}")
    print("\nThe Thurston-Bennequin number is calculated using the formula: tb = w(D) - c_+(D)")
    print(f"tb = {writhe} - {c_plus} = {tb}")
    
solve()