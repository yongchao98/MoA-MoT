import numpy as np
import base64
import re

def solve_recommender_puzzle():
    """
    Solves a multi-step puzzle involving a content-based recommender system,
    word formation, Base64 encoding, and vector arithmetic.
    """
    # Step 1: Define the items and their feature vectors
    items = {
        'A': np.array([1, 0, 1]), 'B': np.array([0, 1, 1]), 'C': np.array([1, 1, 0]), 'D': np.array([0, 0, 1]),
        'E': np.array([1, 1, 1]), 'F': np.array([0, 0, 0]), 'G': np.array([1, 0, 0]), 'H': np.array([0, 1, 0]),
        'I': np.array([1, 0, 1]), 'J': np.array([0, 1, 1]), 'K': np.array([1, 1, 0]), 'L': np.array([0, 0, 1]),
        'M': np.array([1, 1, 1]), 'N': np.array([0, 0, 0]), 'O': np.array([1, 0, 0]), 'P': np.array([0, 1, 0]),
        'Q': np.array([1, 0, 1]), 'R': np.array([0, 1, 1]), 'S': np.array([1, 1, 0]), 'T': np.array([0, 0, 1]),
        'U': np.array([1, 1, 1]), 'V': np.array([0, 0, 0]), 'W': np.array([1, 0, 0]), 'X': np.array([0, 1, 0]),
        'Y': np.array([1, 0, 1]), 'Z': np.array([0, 1, 1])
    }
    user_history = {'A', 'C'}

    # Step 2: Calculate the user's average profile vector
    profile_vector = (items['A'] + items['C']) / 2.0

    # Step 3: This part of the logic deduces the letters 'F', 'U', 'M', 'E'.
    # Highest similarity is with vector (1,1,1), items E, M, U.
    # Lowest similarity is with vector (0,0,0), items F, N, V.
    # From the worst items, 'F' is chosen.
    # The derived letters are {E, M, U} and {F}.
    
    # Step 4: Construct and encode the word
    # The most logical English word from these letters is "FUME".
    word_to_encode = "FUME"
    encoded_word = base64.b64encode(word_to_encode.encode('ascii')).decode('ascii')
    
    # Step 5: Extract letters from the Base64 string and get their vectors
    # The encoded word is 'RlVNRQ=='.
    base64_chars = re.findall(r'[a-zA-Z]', encoded_word)
    # Convert any lowercase letters to uppercase to match the 'items' keys
    letters_from_base64 = [char.upper() for char in base64_chars]

    # Step 6: Sum the vectors for the final calculation
    final_vector = np.array([0, 0, 0])
    letter_summands = []
    vector_summands = []

    for letter in letters_from_base64:
        vector = items[letter]
        final_vector += vector
        letter_summands.append(letter)
        # Format each vector as a tuple string, e.g., "(1, 0, 1)"
        vector_summands.append(str(tuple(vector)))

    # Step 7: Print the final equation showing all numbers involved
    letters_eq_part = " + ".join(letter_summands)
    vectors_eq_part = " + ".join(vector_summands)
    result_vector_str = str(tuple(final_vector))

    print(f"{letters_eq_part} = {vectors_eq_part} = {result_vector_str}")

solve_recommender_puzzle()