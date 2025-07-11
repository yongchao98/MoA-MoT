import numpy as np

def solve_vector_puzzle():
    """
    Solves the vector puzzle by decrypting the message, applying rotations,
    and summing the resulting vectors.
    """
    message = "nggyunglydngraady"

    # 1. Decrypt/filter the message to get a sequence of movements
    direction_map = {
        'n': (0, 1, 0),   # North
        's': (0, -1, 0),  # South
        'e': (1, 0, 0),   # East
        'w': (-1, 0, 0),  # West
        'u': (0, 0, 1),   # Up
        'd': (0, 0, -1),  # Down
        'r': (1, 0, 0),   # Right (East)
        'l': (-1, 0, 0),  # Left (West)
    }

    initial_vectors = []
    for char in message:
        if char in direction_map:
            initial_vectors.append(np.array(direction_map[char]))

    # 2. Apply rotation to every 2nd vector and prepare for summation
    final_vectors = []
    total_sum = np.array([0, 0, 0])

    for i, vec in enumerate(initial_vectors):
        # The problem states "every 2nd vector", which is 1-indexed (2nd, 4th, etc.)
        # This corresponds to 0-indexed list indices 1, 3, 5, etc.
        is_second_vector = (i + 1) % 2 == 0
        
        final_vec = vec
        if is_second_vector:
            # Clockwise rotation around x-axis: (x, y, z) -> (x, z, -y)
            x, y, z = vec
            final_vec = np.array([x, z, -y])
        
        final_vectors.append(final_vec)
        total_sum = total_sum + final_vec

    # 3. Print the full equation and the final sum
    equation_parts = []
    for vec in final_vectors:
        equation_parts.append(f"({vec[0]},{vec[1]},{vec[2]})")
    
    equation_str = " + ".join(equation_parts)
    result_str = f"({total_sum[0]},{total_sum[1]},{total_sum[2]})"
    
    print(f"{equation_str} = {result_str}")

    # Output the final answer in the specified format
    print(f"\n<<<{result_str}>>>")


solve_vector_puzzle()