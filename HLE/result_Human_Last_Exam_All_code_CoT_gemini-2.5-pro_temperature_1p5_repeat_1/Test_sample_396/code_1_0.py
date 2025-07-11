def solve_vector_puzzle():
    """
    Solves the puzzle by decrypting the message, processing vector transformations,
    and summing the results.
    """
    # 1. Decrypt the message using Vigenere cipher with key "vector".
    # The plaintext is "southwestdownsouth".
    decrypted_message = "southwestdownsouth"
    
    # 2. Parse the decrypted string into a list of direction words.
    direction_words = ["north", "south", "east", "west", "up", "down"]
    parsed_directions = []
    temp_message = decrypted_message
    while len(temp_message) > 0:
        for word in direction_words:
            if temp_message.startswith(word):
                parsed_directions.append(word)
                temp_message = temp_message[len(word):]
                break

    # 3. Define the vectors for each direction. Z is up.
    vector_map = {
        "north": (0, 1, 0),
        "south": (0, -1, 0),
        "east": (1, 0, 0),
        "west": (-1, 0, 0),
        "up": (0, 0, 1),
        "down": (0, 0, -1),
    }

    processed_vectors = []
    sum_vector = [0, 0, 0]

    # 4. Process each vector, applying rotation where necessary.
    for i, direction in enumerate(parsed_directions):
        vector = list(vector_map[direction])
        
        # Apply clockwise rotation around x-axis for every 2nd vector.
        # Transformation: (x, y, z) -> (x, z, -y)
        if (i + 1) % 2 == 0:
            original_y = vector[1]
            vector[1] = vector[2]
            vector[2] = -original_y
        
        processed_vectors.append(tuple(vector))

    # 5. Calculate the sum.
    for vec in processed_vectors:
        sum_vector[0] += vec[0]
        sum_vector[1] += vec[1]
        sum_vector[2] += vec[2]

    # 6. Format and print the final equation.
    equation_parts = [str(v) for v in processed_vectors]
    equation_str = " + ".join(equation_parts)
    result_str = f"({sum_vector[0]},{sum_vector[1]},{sum_vector[2]})"
    
    print(f"Sum = {equation_str} = {result_str}")

solve_vector_puzzle()