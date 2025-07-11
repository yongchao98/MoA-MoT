import numpy as np

def solve_and_print():
    """
    Solves the vector decryption and calculation problem.
    """
    message = "nggyunglydngraady"

    # 1. Decrypt the message based on the deduced substitution cipher
    char_to_direction = {
        'n': 'north', 'u': 'up', 'd': 'down', 'l': 'left', 'r': 'right',
        'a': 'east', 'g': 'south', 'y': 'west'
    }
    directions = [char_to_direction[char] for char in message]

    # 2. Define direction to vector mapping
    direction_to_vector = {
        'north': np.array([0, 1, 0]),
        'south': np.array([0, -1, 0]),
        'east': np.array([1, 0, 0]),
        'west': np.array([-1, 0, 0]),
        'up': np.array([0, 0, 1]),
        'down': np.array([0, 0, -1]),
        'left': np.array([-1, 0, 0]),
        'right': np.array([1, 0, 0])
    }

    final_vectors = []
    for i, direction in enumerate(directions):
        vector = direction_to_vector[direction]
        
        # 3. Rotate every 2nd vector (index 1, 3, 5...)
        is_second_vector = (i + 1) % 2 == 0
        if is_second_vector:
            x, y, z = vector
            # Clockwise rotation around x-axis: (x, y, z) -> (x, z, -y)
            vector = np.array([x, z, -y])
            
        final_vectors.append(vector)

    # 4. Calculate the sum and format the output string
    total_sum = np.sum(final_vectors, axis=0)
    
    x_components = [str(v[0]) for v in final_vectors]
    y_components = [str(v[1]) for v in final_vectors]
    z_components = [str(v[2]) for v in final_vectors]

    x_sum_str = "+".join(x_components).replace("+-", "-")
    y_sum_str = "+".join(y_components).replace("+-", "-")
    z_sum_str = "+".join(z_components).replace("+-", "-")

    print(f"({x_sum_str}, {y_sum_str}, {z_sum_str}) = ({total_sum[0]},{total_sum[1]},{total_sum[2]})")

solve_and_print()
<<<(-1,-2,2)>>>