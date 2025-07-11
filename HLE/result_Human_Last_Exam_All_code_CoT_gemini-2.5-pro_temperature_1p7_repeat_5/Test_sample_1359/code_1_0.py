def is_prime(n):
    """
    Checks if a number is prime. For the context of this problem (modulo 12),
    we only need to check for the primes less than 12.
    """
    return n in [2, 3, 5, 7, 11]

def calculate_next_triplet(left_triplet, top_triplet):
    """
    Calculates the next triplet based on the left and top triplets.
    """
    # Unpack triplets
    x_left, y_left, z_left = left_triplet
    x_prev, y_prev, z_prev = top_triplet

    # Step 1: Horizontal Transformation from the left triplet
    # These will be the initial values for the vertical transformation.
    x_interim, y_interim, z_interim = 0, 0, 0
    
    if (x_left + y_left) > 10:
        x_interim = (x_left * 3 - y_left) % 12
        y_interim = (y_left * 2 + 4) % 12
        z_interim = (z_left + x_left) % 12
    else:
        x_interim = (x_left * 2 + y_left) % 12
        y_interim = (y_left * 3 - 2) % 12
        z_interim = (z_left * 2) % 12

    # Step 2: Vertical Transformation using the interim values and the top triplet
    x_final, y_final, z_final = 0, 0, 0
    
    if is_prime(z_prev):
        x_final = (x_interim - 3 + y_prev) % 12
        y_final = (y_interim + x_prev) % 12
        z_final = (z_interim * 2 + x_prev) % 12
    else:
        x_final = (x_interim + 2 - y_prev) % 12
        y_final = (y_interim * 2 - x_prev) % 12
        z_final = (z_interim + 3 + z_prev) % 12
        
    return [x_final, y_final, z_final]

def main():
    # Given known triplets
    M_1_1 = [8, 4, 10]
    M_1_2 = [3, 1, 8]
    M_2_0 = [7, 2, 9]

    # Calculate the first missing triplet M[2][1]
    M_2_1 = calculate_next_triplet(M_2_0, M_1_1)
    print(f"The first missing triplet M[2][1] is: {M_2_1}")

    # Calculate the second missing triplet M[2][2]
    M_2_2 = calculate_next_triplet(M_2_1, M_1_2)
    print(f"The second missing triplet M[2][2] is: {M_2_2}")
    
    # Collect all missing elements
    missing_elements = M_2_1 + M_2_2
    
    # Calculate and print the sum
    total_sum = sum(missing_elements)
    
    # Format the equation string
    equation = " + ".join(map(str, missing_elements))
    
    print("\nCalculating the sum of all missing elements:")
    print(f"{equation} = {total_sum}")

if __name__ == "__main__":
    main()
<<<20>>>