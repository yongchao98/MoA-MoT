import numpy as np

def solve_and_print_equation():
    """
    This script calculates the final vector sum based on the problem description.
    The steps to derive the letters are embedded in the logic:
    1. User profile from A=(1,0,1) and C=(1,1,0) is (1, 0.5, 0.5).
    2. Highest similarity is with E, M, U. Lowest is with F, N, V.
    3. Letters U,M,E and F form the word "FUME".
    4. Base64 encoding of "FUME" is "RlVNRQ==".
    5. This script sums the vectors for the letters R, V, N, R, Q from the encoded string.
    """
    
    # Define the item vectors
    items = {
        'A': np.array([1, 0, 1]), 'B': np.array([0, 1, 1]), 'C': np.array([1, 1, 0]),
        'D': np.array([0, 0, 1]), 'E': np.array([1, 1, 1]), 'F': np.array([0, 0, 0]),
        'G': np.array([1, 0, 0]), 'H': np.array([0, 1, 0]), 'I': np.array([1, 0, 1]),
        'J': np.array([0, 1, 1]), 'K': np.array([1, 1, 0]), 'L': np.array([0, 0, 1]),
        'M': np.array([1, 1, 1]), 'N': np.array([0, 0, 0]), 'O': np.array([1, 0, 0]),
        'P': np.array([0, 1, 0]), 'Q': np.array([1, 0, 1]), 'R': np.array([0, 1, 1]),
        'S': np.array([1, 1, 0]), 'T': np.array([0, 0, 1]), 'U': np.array([1, 1, 1]),
        'V': np.array([0, 0, 0]), 'W': np.array([1, 0, 0]), 'X': np.array([0, 1, 0]),
        'Y': np.array([1, 0, 1]), 'Z': np.array([0, 1, 1])
    }
    
    # The Base64 string is "RlVNRQ=="
    encoded_word = "RlVNRQ=="
    
    # Get the letters from the encoded string that are present in our item set
    letters_to_sum = [char for char in encoded_word if char in items]
    
    # Get the vectors for these letters
    vectors_to_sum = [items[letter] for letter in letters_to_sum]
    
    # Calculate the sum of the vectors
    result_vector = np.sum(vectors_to_sum, axis=0)
    
    # Build the equation string with all numbers
    equation_parts = []
    for vector in vectors_to_sum:
        # Format each vector as (x,y,z) without spaces
        vec_str = str(tuple(vector)).replace(" ", "")
        equation_parts.append(vec_str)
        
    equation_str = " + ".join(equation_parts)
    
    # Format the final result vector
    result_vec_str = str(tuple(result_vector)).replace(" ", "")
    
    # Print the final equation
    print(f"{equation_str} = {result_vec_str}")

# Execute the function
solve_and_print_equation()