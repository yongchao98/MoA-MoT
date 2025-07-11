def is_prime(n):
    """Checks if a number is prime for the purpose of this puzzle."""
    return n in [2, 3, 5, 7, 11]

def calculate_missing_triplets():
    """
    Calculates the missing triplets and the sum of their elements based on the puzzle's rules.
    """
    # Known triplets
    t_3_1 = [7, 2, 9]  # Left of the first missing triplet
    t_2_2 = [8, 4, 10] # Above the first missing triplet
    t_2_3 = [3, 1, 8]  # Above the second missing triplet

    # --- Calculate the first missing triplet: T(3, 2) ---

    # Step 1: Horizontal transformation from T(3, 1)
    xl, yl, zl = t_3_1
    th = [0, 0, 0] # Intermediate triplet
    
    if xl + yl > 10:
        th[0] = (xl * 3 - yl) % 12
        th[1] = (yl * 2 + 4) % 12
        th[2] = (zl + xl) % 12
    else: # xl + yl <= 10
        th[0] = (xl * 2 + yl) % 12
        th[1] = (yl * 3 - 2) % 12
        th[2] = (zl * 2) % 12

    # Step 2: Vertical transformation using Th and T(2, 2)
    xp, yp, zp = t_2_2
    t_3_2 = [0, 0, 0]

    if is_prime(zp):
        t_3_2[0] = (th[0] - 3 + yp) % 12
        t_3_2[1] = (th[1] + xp) % 12
        t_3_2[2] = (th[2] * 2 + xp) % 12
    else: # zp is not prime
        t_3_2[0] = (th[0] + 2 - yp) % 12
        t_3_2[1] = (th[1] * 2 - xp) % 12
        t_3_2[2] = (th[2] + 3 + zp) % 12

    # --- Calculate the second missing triplet: T(3, 3) ---

    # Step 1: Horizontal transformation from the newly found T(3, 2)
    xl, yl, zl = t_3_2
    th2 = [0, 0, 0] # Second intermediate triplet
    
    if xl + yl > 10:
        th2[0] = (xl * 3 - yl) % 12
        th2[1] = (yl * 2 + 4) % 12
        th2[2] = (zl + xl) % 12
    else: # xl + yl <= 10
        th2[0] = (xl * 2 + yl) % 12
        th2[1] = (yl * 3 - 2) % 12
        th2[2] = (zl * 2) % 12
        
    # Step 2: Vertical transformation using Th2 and T(2, 3)
    xp, yp, zp = t_2_3
    t_3_3 = [0, 0, 0]

    if is_prime(zp):
        t_3_3[0] = (th2[0] - 3 + yp) % 12
        t_3_3[1] = (th2[1] + xp) % 12
        t_3_3[2] = (th2[2] * 2 + xp) % 12
    else: # zp is not prime
        t_3_3[0] = (th2[0] + 2 - yp) % 12
        t_3_3[1] = (th2[1] * 2 - xp) % 12
        t_3_3[2] = (th2[2] + 3 + zp) % 12
        
    # --- Sum the elements and print the result ---
    missing_elements = t_3_2 + t_3_3
    total_sum = sum(missing_elements)
    
    equation_parts = [str(e) for e in missing_elements]
    equation = " + ".join(equation_parts)
    
    print(f"The first missing triplet is: {t_3_2}")
    print(f"The second missing triplet is: {t_3_3}")
    print("The sum of the missing elements is calculated as follows:")
    print(f"{equation} = {total_sum}")
    
calculate_missing_triplets()
<<<20>>>