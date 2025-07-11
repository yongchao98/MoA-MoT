import numpy as np

# Helper function for correct modulo operations on negative results
def mod12(n):
    return n % 12

# A simple prime checker for numbers in the range 0-11
def is_prime(n):
    return n in [2, 3, 5, 7, 11]

# Base vertical transformation function based on the puzzle's rules
def V_transform(prev_triplet):
    """
    Applies the base vertical transformation rule to a triplet.
    Assumes Next = f(Previous).
    """
    x_p, y_p, z_p = prev_triplet
    
    if is_prime(z_p):
        # Rule for prime previous z
        next_x = mod12(x_p - 3 + y_p)
        next_y = mod12(y_p + x_p)
        next_z = mod12(z_p * 2 + x_p)
    else:
        # Rule for non-prime previous z
        next_x = mod12(x_p + 2 - y_p)
        next_y = mod12(y_p * 2 - x_p)
        next_z = mod12(z_p + 3 + z_p)
        
    return np.array([next_x, next_y, next_z])

# Given matrix data
T = {
    (1,1): np.array([8, 3, 5]), (1,2): np.array([2, 9, 7]), (1,3): np.array([11, 6, 4]),
    (2,1): np.array([5, 7, 2]), (2,2): np.array([8, 4, 10]),(2,3): np.array([3, 1, 8]),
    (3,1): np.array([7, 2, 9]),
}

# --- Step 1: Calculate the result of V_transform on Row 1 ---
V_res_R2_C1 = V_transform(T[(1,1)])
V_res_R2_C2 = V_transform(T[(1,2)])
V_res_R2_C3 = V_transform(T[(1,3)])

# --- Step 2: Calculate the error vectors for each column ---
# E(c) = Actual_T(2,c) - V_res(2,c)
E1 = mod12(T[(2,1)] - V_res_R2_C1)
E2 = mod12(T[(2,2)] - V_res_R2_C2)
E3 = mod12(T[(2,3)] - V_res_R2_C3)

# --- Step 3: Calculate the missing triplet T(3,2) ---
# It's calculated from T(2,2) with the same logic: V_transform + Error
V_res_R3_C2 = V_transform(T[(2,2)])
T_3_2 = mod12(V_res_R3_C2 + E2)

# --- Step 4: Calculate the missing triplet T(3,3) ---
# It's calculated from T(2,3)
V_res_R3_C3 = V_transform(T[(2,3)])
T_3_3 = mod12(V_res_R3_C3 + E3)

# --- Step 5: Sum the elements and print the result ---
missing_elements = np.concatenate((T_3_2, T_3_3))
total_sum = sum(missing_elements)

# Format the output string as requested
equation_parts = [str(x) for x in missing_elements]
equation_str = " + ".join(equation_parts)

print(f"The first missing triplet (row 3, middle) is: [{T_3_2[0]} {T_3_2[1]} {T_3_2[2]}]")
print(f"The second missing triplet (row 3, right) is: [{T_3_3[0]} {T_3_3[1]} {T_3_3[2]}]")
print("\nThe sum of the missing elements is calculated as follows:")
print(f"{equation_str} = {total_sum}")

print(f"\n<<<The final sum of the missing elements is {total_sum}>>>")
