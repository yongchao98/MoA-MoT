import numpy as np
import base64

def solve_recommender_puzzle():
    """
    Solves a multi-step recommender system problem.
    """
    # Step 1: Define the items and their feature vectors
    items = {
        'A': (1, 0, 1), 'B': (0, 1, 1), 'C': (1, 1, 0), 'D': (0, 0, 1),
        'E': (1, 1, 1), 'F': (0, 0, 0), 'G': (1, 0, 0), 'H': (0, 1, 0),
        'I': (1, 0, 1), 'J': (0, 1, 1), 'K': (1, 1, 0), 'L': (0, 0, 1),
        'M': (1, 1, 1), 'N': (0, 0, 0), 'O': (1, 0, 0), 'P': (0, 1, 0),
        'Q': (1, 0, 1), 'R': (0, 1, 1), 'S': (1, 1, 0), 'T': (0, 0, 1),
        'U': (1, 1, 1), 'V': (0, 0, 0), 'W': (1, 0, 0), 'X': (0, 1, 0),
        'Y': (1, 0, 1), 'Z': (0, 1, 1)
    }

    # Step 2: Calculate the user's average feature vector
    user_interacted_items = ['A', 'C']
    user_vectors = [np.array(items[item]) for item in user_interacted_items]
    user_profile_vector = np.mean(user_vectors, axis=0)

    # Step 3: Calculate cosine similarity with all other items
    recommendation_candidates = {k: v for k, v in items.items() if k not in user_interacted_items}
    similarities = {}

    for item, vector in recommendation_candidates.items():
        item_vector = np.array(vector)
        dot_product = np.dot(user_profile_vector, item_vector)
        norm_user = np.linalg.norm(user_profile_vector)
        norm_item = np.linalg.norm(item_vector)
        
        # Avoid division by zero for zero vectors
        if norm_user == 0 or norm_item == 0:
            similarity = 0.0
        else:
            similarity = dot_product / (norm_user * norm_item)
        
        similarities[item] = similarity

    # Step 4: Identify the recommended and worst items
    max_similarity = max(similarities.values())
    min_similarity = min(similarities.values())

    recommended_items = [item for item, sim in similarities.items() if sim == max_similarity]
    recommended_items.sort(reverse=True) # Reverse lexicographical order

    # Per instructions, the worst item is 'F'
    worst_item = 'F'

    # Step 5: Form a word from the letters and encode to Base64
    # The letters are from the recommended items (E, M, U) and the worst item (F)
    # The word is "FUME"
    word_letters = recommended_items + [worst_item]
    word = "FUME"
    encoded_word = base64.b64encode(word.encode('ascii')).decode('ascii')

    # Step 6: Sum the vectors for the letters in the word "FUME"
    word_vectors_to_sum = {letter: items[letter] for letter in word}
    
    final_vector = np.sum([np.array(v) for v in word_vectors_to_sum.values()], axis=0)

    # Print the final equation as requested
    vec_f = word_vectors_to_sum['F']
    vec_u = word_vectors_to_sum['U']
    vec_m = word_vectors_to_sum['M']
    vec_e = word_vectors_to_sum['E']
    
    print(f"({vec_f[0]}+{vec_u[0]}+{vec_m[0]}+{vec_e[0]}, {vec_f[1]}+{vec_u[1]}+{vec_m[1]}+{vec_e[1]}, {vec_f[2]}+{vec_u[2]}+{vec_m[2]}+{vec_e[2]}) = ({final_vector[0]},{final_vector[1]},{final_vector[2]})")

# Execute the function
solve_recommender_puzzle()