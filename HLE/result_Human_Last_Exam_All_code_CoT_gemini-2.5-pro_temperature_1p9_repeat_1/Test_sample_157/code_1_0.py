import numpy as np
import base64
import operator

def solve_recommendation_task():
    """
    Solves the entire multi-step recommendation problem.
    """
    # Step 1: Define the item vectors
    items = {
        'A': (1, 0, 1), 'B': (0, 1, 1), 'C': (1, 1, 0), 'D': (0, 0, 1), 'E': (1, 1, 1),
        'F': (0, 0, 0), 'G': (1, 0, 0), 'H': (0, 1, 0), 'I': (1, 0, 1), 'J': (0, 1, 1),
        'K': (1, 1, 0), 'L': (0, 0, 1), 'M': (1, 1, 1), 'N': (0, 0, 0), 'O': (1, 0, 0),
        'P': (0, 1, 0), 'Q': (1, 0, 1), 'R': (0, 1, 1), 'S': (1, 1, 0), 'T': (0, 0, 1),
        'U': (1, 1, 1), 'V': (0, 0, 0), 'W': (1, 0, 0), 'X': (0, 1, 0), 'Y': (1, 0, 1),
        'Z': (0, 1, 1)
    }
    item_vectors = {k: np.array(v) for k, v in items.items()}

    # Step 2: Compute the user's profile vector
    user_items = ['A', 'C']
    vector_A = item_vectors['A']
    vector_C = item_vectors['C']
    user_profile_vector = (vector_A + vector_C) / 2
    print(f"The user has interacted with A={tuple(vector_A)} and C={tuple(vector_C)}.")
    print(f"The user's average profile vector is ({vector_A[0]}+{vector_C[0]})/2, ({vector_A[1]}+{vector_C[1]})/2, ({vector_A[2]}+{vector_C[2]})/2 = {tuple(user_profile_vector)}\n")

    # Step 3: Calculate cosine similarity for all other items
    def cosine_similarity(v1, v2):
        norm_v1 = np.linalg.norm(v1)
        norm_v2 = np.linalg.norm(v2)
        if norm_v1 == 0 or norm_v2 == 0:
            return 0.0
        return np.dot(v1, v2) / (norm_v1 * norm_v2)

    recommendation_candidates = {k: v for k, v in item_vectors.items() if k not in user_items}
    similarities = {k: cosine_similarity(user_profile_vector, v) for k, v in recommendation_candidates.items()}

    # Step 4: Find the recommended and worst items
    max_similarity = max(similarities.values())
    min_similarity = min(similarities.values())

    recommended_items = sorted([k for k, v in similarities.items() if np.isclose(v, max_similarity)], reverse=True)
    worst_item = 'F' # "the one students in the USA would least like to get on a test"

    print(f"The item(s) with the highest similarity are {', '.join(recommended_items)}.")
    print(f"The recommended items in reverse lexicographical order are: {','.join(recommended_items)}")
    print(f"The worst item (lowest similarity and hint) is: {worst_item}\n")

    # Step 5: Construct and encode the word
    # The letters from the previous step are U, M, E and F. An English word from these is FUME.
    word_to_encode = "FUME"
    encoded_word = base64.b64encode(word_to_encode.encode('ascii')).decode('ascii')
    print(f"Using the letters {','.join(recommended_items)} and {worst_item} to form the word: {word_to_encode}")
    print(f"Encoding '{word_to_encode}' in Base64 results in: {encoded_word}\n")

    # Step 6: Extract letters from Base64 string and handle cases
    base64_raw_letters = [char for char in encoded_word if char.isalpha()]
    # All items are uppercase, so we convert the lowercase 'l' to 'L'
    final_letters = [letter.upper() for letter in base64_raw_letters]
    print(f"The letters extracted from the Base64 string are: {', '.join(final_letters)}\n")
    
    # Step 7: Sum the corresponding vectors and show the equation
    vectors_to_sum = [item_vectors[letter] for letter in final_letters]
    
    x_components = [str(v[0]) for v in vectors_to_sum]
    y_components = [str(v[1]) for v in vectors_to_sum]
    z_components = [str(v[2]) for v in vectors_to_sum]
    
    final_vector = sum(vectors_to_sum)
    
    print("Adding the vectors for these letters:")
    vector_strings = [str(tuple(v)) for v in vectors_to_sum]
    print(f"  {' + '.join(vector_strings)}")
    print("Calculating each component of the resulting vector:")
    print(f"  x = {' + '.join(x_components)} = {final_vector[0]}")
    print(f"  y = {' + '.join(y_components)} = {final_vector[1]}")
    print(f"  z = {' + '.join(z_components)} = {final_vector[2]}\n")

    # Final Answer
    final_answer_tuple = tuple(final_vector)
    print(f"The final resulting vector is {final_answer_tuple}")
    print(f"<<<{final_answer_tuple}>>>")

solve_recommendation_task()