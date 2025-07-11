def is_prime(n):
    """Checks if a number is prime (for numbers up to 11)."""
    return n in [2, 3, 5, 7, 11]

def horizontal_transform(L):
    """Applies horizontal transformation rules to a triplet L."""
    x, y, z = L
    if x + y > 10:
        h_x = (x * 3 - y) % 12
        h_y = (y * 2 + 4) % 12
        h_z = (z + x) % 12
    else:
        h_x = (x * 2 + y) % 12
        h_y = (y * 3 - 2) % 12
        h_z = (z * 2) % 12
    return [h_x, h_y, h_z]

def vertical_self_transform(U):
    """Applies vertical transformation rules to a triplet U, using U as its own 'previous'."""
    x, y, z = U
    prev_x, prev_y, prev_z = U
    if is_prime(prev_z):
        v_x = (x - 3 + prev_y) % 12
        v_y = (y + prev_x) % 12
        v_z = (z * 2 + prev_x) % 12
    else:
        v_x = (x + 2 - prev_y) % 12
        v_y = (y * 2 - prev_x) % 12
        v_z = (z + 3 + prev_z) % 12
    return [v_x, v_y, v_z]

def solve_missing_triplets():
    """
    Calculates the missing triplets based on the derived complex interaction model.
    """
    T31 = [7, 2, 9]
    T22 = [8, 4, 10]
    T23 = [3, 1, 8]

    # --- Calculate T32 ---
    L = T31
    U = T22

    # Calculate x32 (middle column logic)
    H_L_for_32 = horizontal_transform(L)
    x32 = H_L_for_32[0]

    # Calculate y32 and z32
    V_U_for_32 = vertical_self_transform(U)
    if is_prime(U[2]):
        # This branch is not taken for T32, as U[2]=10 is not prime
        y32 = (V_U_for_32[1] - L[1]) % 12
        z32 = (V_U_for_32[2] + L[0] + 1) % 12
    else: # U[2]=10 is not prime
        y32 = V_U_for_32[1]
        z32 = (V_U_for_32[2] - 3) % 12
        
    T32 = [x32, y32, z32]
    
    # --- Calculate T33 ---
    L = T32
    U = T23

    # Calculate x33 (right column logic)
    H_L_for_33 = horizontal_transform(L)
    y_middle = T32[1] # y_3,2
    x33 = (H_L_for_33[0] - y_middle - 1) % 12

    # Calculate y33 and z33
    V_U_for_33 = vertical_self_transform(U)
    if is_prime(U[2]):
         # This branch is not taken for T33, as U[2]=8 is not prime
        y33 = (V_U_for_33[1] - L[1]) % 12
        z33 = (V_U_for_33[2] + L[0] + 1) % 12
    else: # U[2]=8 is not prime
        y33 = V_U_for_33[1]
        z33 = (V_U_for_33[2] - 3) % 12

    T33 = [x33, y33, z33]

    # --- Print results and sum ---
    print(f"The missing triplet at position (3, 2) is: {T32}")
    print(f"The missing triplet at position (3, 3) is: {T33}")
    
    total_sum = sum(T32) + sum(T33)
    
    print("\nThe sum of the missing elements is calculated as follows:")
    sum_equation = " + ".join(map(str, T32 + T33))
    print(f"{sum_equation} = {total_sum}")

solve_missing_triplets()
<<<34>>>