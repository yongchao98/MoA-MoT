def is_prime(n):
    """Checks if a number is prime for the context of this problem."""
    return n in [2, 3, 5, 7, 11]

def horizontal_transform(triplet):
    """Applies the horizontal transformation rules."""
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

def vertical_transform(prev_triplet):
    """Applies the vertical transformation rules."""
    x, y, z = prev_triplet
    if is_prime(z):
        # Rule: Next x = (x - 3 + previous y) mod 12. Here 'x' is prev_x.
        next_x = (x - 3 + y) % 12
        next_y = (y + x) % 12
        next_z = (z * 2 + x) % 12
    else:
        # Rule: Next x = (x + 2 - previous y) mod 12. Here 'x' is prev_x.
        next_x = (x + 2 - y) % 12
        next_y = (y * 2 - x) % 12
        next_z = (z + 3 + z) % 12
    return [next_x, next_y, next_z]

def calculate_next_triplet(left_triplet, top_triplet):
    """Calculates a triplet based on its left and top neighbors."""
    h_trans = horizontal_transform(left_triplet)
    v_trans = vertical_transform(top_triplet)

    # Combination logic derived from analyzing the grid
    final_x = h_trans[0]
    final_y = (h_trans[1] + v_trans[1] - 1) % 12
    final_z = (h_trans[2] + v_trans[2] - 1) % 12
    
    return [final_x, final_y, final_z]

def solve_puzzle():
    """Main function to solve the puzzle."""
    # Known values from the grid
    l3 = [7, 2, 9]
    m2 = [8, 4, 10]
    r2 = [3, 1, 8]

    # Calculate the missing middle triplet in the last row
    m3 = calculate_next_triplet(l3, m2)
    
    # Calculate the missing right triplet in the last row
    r3 = calculate_next_triplet(m3, r2)
    
    missing_numbers = m3 + r3
    total_sum = sum(missing_numbers)
    
    equation_str = " + ".join(map(str, missing_numbers))
    
    print(f"The first missing triplet is: {m3}")
    print(f"The second missing triplet is: {r3}")
    print(f"The sum of the missing elements is calculated as follows:")
    print(f"{equation_str} = {total_sum}")

solve_puzzle()
<<<29>>>