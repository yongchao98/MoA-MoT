import math

def is_prime(n):
    """Checks if a number is prime (for the context of this problem)."""
    return n in [2, 3, 5, 7, 11]

def horizontal_transform(triplet):
    """Applies the horizontal transformation rule."""
    x, y, z = triplet
    if x + y > 10:
        next_x = (x * 3 - y) % 12
        next_y = (y * 2 + 4) % 12
        next_z = (z + x) % 12
    else:
        next_x = (x * 2 + y) % 12
        next_y = (y * 3 - 2) % 12
        next_z = (z * 2) % 12
    return [next_x, next_y, next_z]

def vertical_transform(base_triplet, prev_triplet):
    """Applies the vertical transformation rule."""
    x, y, z = base_triplet
    prev_x, prev_y, prev_z = prev_triplet

    if is_prime(prev_z):
        next_x = (x - 3 + prev_y) % 12
        next_y = (y + prev_x) % 12
        next_z = (z * 2 + prev_x) % 12
    else:
        next_x = (x + 2 - prev_y) % 12
        next_y = (y * 2 - prev_x) % 12
        next_z = (z + 3 + prev_z) % 12
    return [next_x, next_y, next_z]

def solve_puzzle():
    """Main function to solve the puzzle."""
    
    # Given triplets from the matrix
    t31 = [7, 2, 9]
    t22 = [8, 4, 10]
    t23 = [3, 1, 8]
    
    # --- Calculate the first missing triplet T(3,2) ---
    # 1. Horizontal transformation of T(3,1)
    h_out_1 = horizontal_transform(t31)
    # 2. Vertical transformation using the result and T(2,2)
    t32 = vertical_transform(h_out_1, t22)
    
    # --- Calculate the second missing triplet T(3,3) ---
    # 1. Horizontal transformation of the newly found T(3,2)
    h_out_2 = horizontal_transform(t32)
    # 2. Vertical transformation using the result and T(2,3)
    t33_pre_adjustment = vertical_transform(h_out_2, t23)
    
    # 3. Apply the cross-dependency rule for the rightmost column
    t33 = t33_pre_adjustment
    # The middle triplet is T(3,2), its y value influences T(3,3)'s x value.
    y_middle = t32[1] 
    adjustment = round(y_middle / 4.0)
    t33[0] = (t33[0] - adjustment + 12) % 12 # add 12 to handle potential negatives
    
    # The missing elements are the components of t32 and t33
    missing_elements = t32 + t33
    total_sum = sum(missing_elements)
    
    # Print the final equation
    equation_str = " + ".join(map(str, missing_elements))
    print(f"The missing triplet [? ? ?] in the middle of the last row is {t32}.")
    print(f"The missing triplet [? ? ?] on the right of the last row is {t33}.")
    print("\nThe sum of the missing elements is calculated as follows:")
    print(f"{equation_str} = {total_sum}")
    
    return total_sum

# Run the solver and output the result in the specified format.
final_sum = solve_puzzle()
# The final answer tag will be generated after the print statements.
# No extra code here, the helper function takes care of the output.
print(f'<<<{final_sum}>>>')