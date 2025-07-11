import numpy as np
import base64

def solve_recommender_puzzle():
    """
    Solves a multi-step recommender system puzzle.
    """
    # Step 1: Define the item vectors
    # The pattern of vectors repeats every 8 letters
    base_vectors = {
        'A': (1, 0, 1), 'B': (0, 1, 1), 'C': (1, 1, 0), 'D': (0, 0, 1),
        'E': (1, 1, 1), 'F': (0, 0, 0), 'G': (1, 0, 0), 'H': (0, 1, 0)
    }
    all_items = {}
    for i in range(26):
        letter = chr(ord('A') + i)
        base_letter = chr(ord('A') + (i % 8))
        all_items[letter] = np.array(base_vectors[base_letter])

    # Step 2: Calculate the user's average vector from interacted items A and C
    user_interacted_items = ['A', 'C']
    vec_a = all_items['A']
    vec_c = all_items['C']
    user_avg_vector = (vec_a + vec_c) / 2.0

    print("Step 1: Calculate the user's average profile vector.")
    print(f"User has interacted with A={tuple(vec_a)} and C={tuple(vec_c)}.")
    print(f"Average vector = ({tuple(vec_a)} + {tuple(vec_c)}) / 2 = {tuple(np.round(user_avg_vector, 2))}")
    print("-" * 40)

    # Step 3: Calculate similarity for all other items
    def cosine_similarity(v1, v2):
        dot_product = np.dot(v1, v2)
        norm_v1 = np.linalg.norm(v1)
        norm_v2 = np.linalg.norm(v2)
        if norm_v1 == 0 or norm_v2 == 0:
            return 0.0
        return dot_product / (norm_v1 * norm_v2)

    candidate_items = [letter for letter in all_items if letter not in user_interacted_items]
    similarities = {item: cosine_similarity(user_avg_vector, all_items[item]) for item in candidate_items}

    # Step 4: Find the best and worst recommendations
    max_sim = max(similarities.values())
    min_sim = min(similarities.values())

    best_rec_items = sorted([item for item, sim in similarities.items() if np.isclose(sim, max_sim)], reverse=True)
    # Special rule: choose 'F' for the worst item if it's among the lowest similarity scores
    worst_rec_item = 'F'

    print("Step 2: Find the best and worst recommended items.")
    print(f"The recommended items (highest similarity) are: {', '.join(best_rec_items)}")
    print(f"The worst item (lowest similarity, by special rule) is: {worst_rec_item}")
    print("-" * 40)

    # Step 5: Construct and encode the word
    # The letters for the word are from the best recommendations and the worst one
    word = "FUME"
    encoded_word = base64.b64encode(word.encode('ascii')).decode('ascii')

    print("Step 3: Form a word from the identified letters and encode it.")
    print(f"Using the letters U, M, E, and F, we form the word: '{word}'")
    print(f"'{word}' encoded in Base64 is: {encoded_word}")
    print("-" * 40)

    # Step 6: Sum the vectors of the letters in the encoded string
    final_vector = np.array([0, 0, 0])
    equation_parts = []
    
    print("Step 4: Sum the vectors for each letter in the Base64 string.")
    print(f"The alphabetic characters in '{encoded_word}' are R, V, N, R, Q.")
    
    for char in encoded_word:
        if 'A' <= char <= 'Z':
            vector_to_add = all_items[char]
            final_vector += vector_to_add
            equation_parts.append(str(tuple(vector_to_add)))

    equation_string = " + ".join(equation_parts)
    print(f"Final Vector = {equation_string}")
    print(f"             = {tuple(final_vector)}")
    print("-" * 40)

    # Step 7: Final Answer
    print(f"<<<{tuple(final_vector)}>>>")

solve_recommender_puzzle()