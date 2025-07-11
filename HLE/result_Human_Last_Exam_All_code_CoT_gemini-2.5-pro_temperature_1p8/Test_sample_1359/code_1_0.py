def solve_puzzle():
    """
    Calculates the missing elements in a 3x3 matrix based on a set of rules
    and then finds their sum.
    """

    # Prime numbers in the range 0-11
    primes = {2, 3, 5, 7, 11}

    # Known matrix values
    m20 = [7, 2, 9]  # M[2][0]
    m11 = [8, 4, 10] # M[1][1]
    m12 = [3, 1, 8]  # M[1][2]

    # --- Step 1: Calculate M[2][1] ---

    # A. Horizontal transformation from M[2][0] = [7, 2, 9]
    h_prev_x, h_prev_y, h_prev_z = m20
    
    # Check condition for horizontal rule
    if h_prev_x + h_prev_y > 10:
        # This branch is not taken for M[2][0]
        pass
    else: # 7 + 2 <= 10
        interim_x = (h_prev_x * 2 + h_prev_y) % 12
        interim_y = (h_prev_y * 3 - 2) % 12
        interim_z = (h_prev_z * 2) % 12
    
    # Intermediate triplet is [16%12, 4%12, 18%12] = [4, 4, 6]

    # B. Vertical transformation using M[1][1] = [8, 4, 10]
    v_prev_x, v_prev_y, v_prev_z = m11
    
    # Check condition for vertical rule
    if v_prev_z in primes:
        # This branch is not taken for M[1][1]
        pass
    else: # 10 is not prime
        final_x1 = (interim_x + 2 - v_prev_y) % 12
        final_y1 = (interim_y * 2 - v_prev_x) % 12
        final_z1 = (interim_z + 3 + v_prev_z) % 12

    m21 = [final_x1, final_y1, final_z1]
    # m21 is [(4+2-4)%12, (4*2-8)%12, (6+3+10)%12] = [2, 0, 7]

    # --- Step 2: Calculate M[2][2] ---
    
    # A. Horizontal transformation from the new M[2][1] = [2, 0, 7]
    h_prev_x, h_prev_y, h_prev_z = m21
    
    # Check condition for horizontal rule
    if h_prev_x + h_prev_y > 10:
        # This branch is not taken for M[2][1]
        pass
    else: # 2 + 0 <= 10
        interim_x = (h_prev_x * 2 + h_prev_y) % 12
        interim_y = (h_prev_y * 3 - 2) % 12
        interim_z = (h_prev_z * 2) % 12

    # Intermediate triplet is [4%12, -2%12, 14%12] = [4, 10, 2]

    # B. Vertical transformation using M[1][2] = [3, 1, 8]
    v_prev_x, v_prev_y, v_prev_z = m12
    
    # Check condition for vertical rule
    if v_prev_z in primes:
         # This branch is not taken for M[1][2]
        pass
    else: # 8 is not prime
        final_x2 = (interim_x + 2 - v_prev_y) % 12
        final_y2 = (interim_y * 2 - v_prev_x) % 12
        final_z2 = (interim_z + 3 + v_prev_z) % 12

    m22 = [final_x2, final_y2, final_z2]
    # m22 is [(4+2-1)%12, (10*2-3)%12, (2+3+8)%12] = [5, 5, 1]
    
    # --- Step 3: Calculate and print the sum ---
    
    total_sum = sum(m21) + sum(m22)

    # Print the equation as requested
    print(f"{m21[0]} + {m21[1]} + {m21[2]} + {m22[0]} + {m22[1]} + {m22[2]} = {total_sum}")

solve_puzzle()
<<<20>>>