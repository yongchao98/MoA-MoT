import sys

def is_prime(n):
    """Checks if a number is prime for modulo 12 context."""
    return n in [2, 3, 5, 7, 11]

def calculate_triplets():
    """
    Calculates the missing triplets based on the specified rules and prints the results.
    """
    # Known triplets
    T_31 = [7, 2, 9]
    T_22 = [8, 4, 10]
    T_23 = [3, 1, 8]
    
    missing_elements = []
    
    # --- Calculate T_32 ---
    print("Calculating Triplet at Row 3, Column 2:")
    
    # Step 1: Horizontal transformation from T_left = T_31
    t_left = T_31
    x, y, z = t_left
    
    # x+y = 7+2=9 <= 10
    h_x = (x * 2 + y) % 12
    h_y = (y * 3 - 2) % 12
    h_z = (z * 2) % 12
    t_int = [h_x, h_y, h_z]
    
    # Step 2: Vertical transformation with t_int and T_top = T_22
    t_curr = t_int
    t_prev = T_22
    x_c, y_c, z_c = t_curr
    x_p, y_p, z_p = t_prev

    # previous z = 10 (not prime)
    x_new = (x_c + 2 - y_p) % 12
    y_new = (y_c * 2 - x_p) % 12
    z_new = (z_c + 3 + z_p) % 12
    
    T_32 = [x_new, y_new, z_new]
    missing_elements.extend(T_32)
    
    print(f"H-Step on {t_left} gave intermediate {t_int}")
    print(f"V-Step on {t_int} with top {t_prev} gave final {T_32}")
    print(f"x = ({x_c} + 2 - {y_p}) % 12 = {T_32[0]}")
    print(f"y = ({y_c} * 2 - {x_p}) % 12 = {T_32[1]}")
    print(f"z = ({z_c} + 3 + {z_p}) % 12 = {T_32[2]}\n")

    # --- Calculate T_33 ---
    print("Calculating Triplet at Row 3, Column 3:")
    
    # Step 1: Horizontal transformation from T_left = T_32
    t_left = T_32
    x, y, z = t_left

    # x+y = T_32[0]+T_32[1] = 2+0=2 <= 10
    h_x = (x * 2 + y) % 12
    h_y = (y * 3 - 2) % 12
    h_z = (z * 2) % 12
    t_int = [h_x, h_y, h_z]
    
    # Step 2: Vertical transformation with t_int and T_top = T_23
    t_curr = t_int
    t_prev = T_23
    x_c, y_c, z_c = t_curr
    x_p, y_p, z_p = t_prev

    # previous z = 8 (not prime)
    x_new = (x_c + 2 - y_p) % 12
    y_new = (y_c * 2 - x_p) % 12
    z_new = (z_c + 3 + z_p) % 12
    
    T_33 = [x_new, y_new, z_new]
    missing_elements.extend(T_33)

    print(f"H-Step on {t_left} gave intermediate {t_int}")
    print(f"V-Step on {t_int} with top {t_prev} gave final {T_33}")
    print(f"x = ({x_c} + 2 - {y_p}) % 12 = {T_33[0]}")
    print(f"y = ({y_c} * 2 - {x_p}) % 12 = {T_33[1]}")
    print(f"z = ({z_c} + 3 + {z_p}) % 12 = {T_33[2]}\n")
    
    # --- Summation ---
    total_sum = sum(missing_elements)
    print("The missing triplets are:", T_32, "and", T_33)
    print("The sum of their elements is:", total_sum)

calculate_triplets()