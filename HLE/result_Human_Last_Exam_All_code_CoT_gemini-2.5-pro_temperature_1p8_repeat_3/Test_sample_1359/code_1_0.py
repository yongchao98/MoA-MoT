def is_prime(n):
    """Checks if a number is prime for this puzzle's context."""
    return n in [2, 3, 5, 7, 11]

def vertical_transform(triplet):
    """Applies the vertical transformation rules to a triplet."""
    x, y, z = triplet
    if is_prime(z):
        # If previous z is prime
        next_x = (x - 3 + y) % 12
        next_y = (y + x) % 12
        next_z = (z * 2 + x) % 12
    else:
        # If previous z is not prime
        next_x = (x + 2 - y) % 12
        next_y = (y * 2 - x) % 12
        next_z = (z * 2 + 3) % 12
    return [next_x, next_y, next_z]

def horizontal_transform(triplet):
    """Applies the horizontal transformation rules to a triplet."""
    x, y, z = triplet
    if x + y > 10:
        # If x + y > 10
        next_x = (x * 3 - y) % 12
        next_y = (y * 2 + 4) % 12
        next_z = (z + x) % 12
    else:
        # If x + y <= 10
        next_x = (x * 2 + y) % 12
        next_y = (y * 3 - 2) % 12
        next_z = (z * 2) % 12
    return [next_x, next_y, next_z]

# Given triplets
T_2_1 = [5, 7, 2]
T_2_2 = [8, 4, 10]
T_2_3 = [3, 1, 8]
T_3_1 = [7, 2, 9]

# --- Calculate T(3,2) ---
# V from T(2,2)
V_from_T22 = vertical_transform(T_2_2)
# H from T(3,1)
H_from_T31 = horizontal_transform(T_3_1)

# Combination rule for column j=2
x_32 = V_from_T22[0] # v_x
y_32 = V_from_T22[2] # v_z
z_32 = (H_from_T31[1] + V_from_T22[2]) % 12 # h_y + v_z

T_3_2 = [x_32, y_32, z_32]

# --- Calculate T(3,3) ---
# V from T(2,3)
V_from_T23 = vertical_transform(T_2_3)
# H from T(3,2)
H_from_T32 = horizontal_transform(T_3_2)

# Combination rule for column j=3
x_33 = (V_from_T23[2] - H_from_T32[0]) % 12 # v_z - h_x
y_33 = V_from_T23[1] # v_y
z_33 = H_from_T32[0] # h_x

T_3_3 = [x_33, y_33, z_33]

# --- Calculate and print the sum ---
missing_elements = T_3_2 + T_3_3
total_sum = sum(missing_elements)

# Print the equation as requested
equation_str = " + ".join(map(str, missing_elements))
print(f"The missing elements for the two triplets are {T_3_2} and {T_3_3}.")
print(f"The sum is calculated as: {equation_str} = {total_sum}")

print(f"\nFinal Answer: The sum of the missing elements is {total_sum}")

<<<38>>>