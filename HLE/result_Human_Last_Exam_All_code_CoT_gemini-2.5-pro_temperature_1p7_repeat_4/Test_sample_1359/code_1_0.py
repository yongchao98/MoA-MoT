import sys

# On some platforms, negative modulo is handled differently.
# This ensures (-a) % n gives a positive result.
def proper_modulo(a, n):
    return a % n

def solve():
    """
    This function calculates the missing elements and their sum based on the puzzle's rules.
    """
    # Given known triplets
    T31 = [7, 2, 9]
    T22 = [8, 4, 10]
    T23 = [3, 1, 8]

    # --- Step 1: Calculate T(3, 2) ---

    # 1a: Horizontal transform from T(3, 1)
    x, y, z = T31
    # x + y = 9 <= 10, so use H-Rule 2
    h_next_x = proper_modulo((x * 2 + y), 12)
    h_next_y = proper_modulo((y * 3 - 2), 12)
    h_next_z = proper_modulo((z * 2), 12)
    
    intermediate_T32 = [h_next_x, h_next_y, h_next_z]

    # 1b: Vertical transform using T(2, 2)
    x, y, z = intermediate_T32
    prev_x, prev_y, prev_z = T22
    
    # prev_z = 10 is not prime, use V-Rule 2
    T32_x = proper_modulo((x + 2 - prev_y), 12)
    T32_y = proper_modulo((y * 2 - prev_x), 12)
    T32_z = proper_modulo((z + 3 + prev_z), 12)
    
    T32 = [T32_x, T32_y, T32_z]

    # --- Step 2: Calculate T(3, 3) ---

    # 2a: Horizontal transform from T(3, 2)
    x, y, z = T32
    # x + y = 2 + 0 = 2 <= 10, so use H-Rule 2
    h_next_x = proper_modulo((x * 2 + y), 12)
    h_next_y = proper_modulo((y * 3 - 2), 12)
    h_next_z = proper_modulo((z * 2), 12)
    
    intermediate_T33 = [h_next_x, h_next_y, h_next_z]

    # 2b: Vertical transform using T(2, 3)
    x, y, z = intermediate_T33
    prev_x, prev_y, prev_z = T23

    # prev_z = 8 is not prime, use V-Rule 2
    T33_x = proper_modulo((x + 2 - prev_y), 12)
    T33_y = proper_modulo((y * 2 - prev_x), 12)
    T33_z = proper_modulo((z + 3 + prev_z), 12)

    T33 = [T33_x, T33_y, T33_z]

    # --- Step 3: Calculate the sum ---
    missing_elements = T32 + T33
    total_sum = sum(missing_elements)
    
    print(f"The first missing triplet is [{T32[0]}, {T32[1]}, {T32[2]}]")
    print(f"The second missing triplet is [{T33[0]}, {T33[1]}, {T33[2]}]")
    equation = " + ".join(map(str, missing_elements))
    print(f"The sum is calculated as: {equation} = {total_sum}")

solve()
<<<20>>>