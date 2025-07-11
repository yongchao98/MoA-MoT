import collections

def solve_vector_puzzle():
    """
    Solves the vector puzzle by decrypting the message, processing the directions,
    and calculating the sum of the resulting vectors.
    """
    # 1. Decrypt the message
    ciphertext = "nggyunglydngraady"
    shift = -2
    
    decrypted_message = ""
    for char in ciphertext:
        decrypted_char_code = (ord(char) - ord('a') + shift) % 26
        decrypted_message += chr(decrypted_char_code + ord('a'))
    
    # 2. Extract direction components
    direction_map = {
        'n': "North", 's': "South", 'e': "East", 'w': "West", 'u': "Up", 'd': "Down"
    }
    
    directions = []
    for char in decrypted_message:
        if char in direction_map:
            directions.append(char)
            
    # 3. Convert directions to vectors
    vector_map = {
        'e': (1, 0, 0),
        'w': (-1, 0, 0),
        'n': (0, 1, 0),
        's': (0, -1, 0),
        'u': (0, 0, 1),
        'd': (0, 0, -1),
    }
    
    initial_vectors = [vector_map[d] for d in directions]
    
    # 4. Apply rotation to every 2nd vector
    final_vectors = []
    for i, vec in enumerate(initial_vectors):
        # Index i is 0-based, so we check i+1 for the "2nd, 4th, ..." vectors
        if (i + 1) % 2 == 0:
            x, y, z = vec
            # Clockwise rotation around x-axis: (x, y, z) -> (x, z, -y)
            rotated_vec = (x, z, -y)
            final_vectors.append(rotated_vec)
        else:
            final_vectors.append(vec)
            
    # 5. Sum the final vectors
    sum_x = sum(v[0] for v in final_vectors)
    sum_y = sum(v[1] for v in final_vectors)
    sum_z = sum(v[2] for v in final_vectors)
    
    # 6. Output the result, showing the equation for each component
    x_components = [str(v[0]) for v in final_vectors]
    y_components = [str(v[1]) for v in final_vectors]
    z_components = [str(v[2]) for v in final_vectors]

    print("Final vector sum calculation:")
    print(f"x = {' + '.join(x_components).replace('+ -', '- ')} = {sum_x}")
    print(f"y = {' + '.join(y_components).replace('+ -', '- ')} = {sum_y}")
    print(f"z = {' + '.join(z_components).replace('+ -', '- ')} = {sum_z}")

    final_sum_vector = (sum_x, sum_y, sum_z)
    print(f"\nThe sum of all generated vectors is {final_sum_vector}")


solve_vector_puzzle()
print("<<<(1, 0, 1)>>>")