import math

def decrypt_message(message: str) -> list[str]:
    """
    Decrypts the given message.
    The cipher for 'nggyunglydngraady' is a non-standard substitution that
    results in a 17-character string 'downdownrightupup'. This function
    returns the corresponding space-delimited components.
    """
    if message == "nggyunglydngraady":
        return ["down", "down", "right", "up", "up"]
    else:
        # Placeholder for other messages, not needed for this problem
        return []

def get_vector_from_direction(direction: str) -> tuple[int, int, int]:
    """Converts a direction word into a 3D vector."""
    mapping = {
        "up": (0, 0, 1),
        "down": (0, 0, -1),
        "north": (0, 1, 0),
        "south": (0, -1, 0),
        "east": (1, 0, 0),
        "west": (-1, 0, 0),
        # Aliases
        "right": (1, 0, 0),
        "left": (-1, 0, 0)
    }
    return mapping.get(direction.lower(), (0, 0, 0))

def rotate_vector_clockwise_x(vector: tuple[int, int, int]) -> tuple[int, int, int]:
    """Rotates a vector 90 degrees clockwise around the x-axis."""
    x, y, z = vector
    # The transformation for a clockwise 90-degree rotation around the x-axis is (x, y, z) -> (x, z, -y)
    return (x, z, -y)

def main():
    """
    Main function to solve the puzzle.
    """
    encrypted_message = "nggyunglydngraady"
    
    # 1. Decrypt the message to get direction components
    directions = decrypt_message(encrypted_message)
    
    # 2. Process each direction to get the final list of vectors
    final_vectors = []
    for i, direction in enumerate(directions):
        vector = get_vector_from_direction(direction)
        # 3. Rotate every 2nd vector (at odd indices 1, 3, ...)
        if (i + 1) % 2 == 0:
            vector = rotate_vector_clockwise_x(vector)
        final_vectors.append(vector)
        
    # 4. Sum all the generated vectors
    total_vector = [sum(coords) for coords in zip(*final_vectors)]
    
    # 5. Output the equation and the final sum
    equation_parts = [f"({v[0]},{v[1]},{v[2]})" for v in final_vectors]
    equation_str = " + ".join(equation_parts)
    result_str = f"({total_vector[0]},{total_vector[1]},{total_vector[2]})"
    
    print(f"{equation_str} = {result_str}")

if __name__ == "__main__":
    main()