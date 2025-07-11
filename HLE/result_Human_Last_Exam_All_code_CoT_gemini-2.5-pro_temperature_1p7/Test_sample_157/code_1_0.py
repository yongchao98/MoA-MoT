import numpy as np
import base64
import math

def solve_recommender_puzzle():
    """
    Solves a multi-step recommender system problem.
    """
    items = {
        'A': (1, 0, 1), 'B': (0, 1, 1), 'C': (1, 1, 0), 'D': (0, 0, 1),
        'E': (1, 1, 1), 'F': (0, 0, 0), 'G': (1, 0, 0), 'H': (0, 1, 0),
        'I': (1, 0, 1), 'J': (0, 1, 1), 'K': (1, 1, 0), 'L': (0, 0, 1),
        'M': (1, 1, 1), 'N': (0, 0, 0), 'O': (1, 0, 0), 'P': (0, 1, 0),
        'Q': (1, 0, 1), 'R': (0, 1, 1), 'S': (1, 1, 0), 'T': (0, 0, 1),
        'U': (1, 1, 1), 'V': (0, 0, 0), 'W': (1, 0, 0), 'X': (0, 1, 0),
        'Y': (1, 0, 1), 'Z': (0, 1, 1)
    }

    # Step 1: Calculate user profile vector
    user_items = ['A', 'C']
    user_vectors = [np.array(items[item]) for item in user_items]
    user_profile = np.mean(user_vectors, axis=0)

    # Step 2: Calculate cosine similarities for candidate items
    similarities = {}
    candidate_items = {k: v for k, v in items.items() if k not in user_items}

    for item, vector in candidate_items.items():
        vec_np = np.array(vector)
        dot_product = np.dot(user_profile, vec_np)
        norm_user = np.linalg.norm(user_profile)
        norm_item = np.linalg.norm(vec_np)

        if norm_user == 0 or norm_item == 0:
            similarity = 0.0
        else:
            similarity = dot_product / (norm_user * norm_item)
        
        similarities[item] = similarity

    # Step 3: Find best and worst items
    max_similarity = max(similarities.values())
    best_items = sorted([item for item, sim in similarities.items() if math.isclose(sim, max_similarity)], reverse=True)

    min_similarity = min(similarities.values())
    worst_items = [item for item, sim in similarities.items() if math.isclose(sim, min_similarity)]
    worst_item = 'F' if 'F' in worst_items else sorted(worst_items)[0]

    recommended_str = ",".join(best_items)
    print(f"{recommended_str} {worst_item}")
    
    # Step 4: Construct word, encode, and sum vectors
    # The letters from the best (U, M, E) and worst (F) items form the word "FUME"
    word_to_encode = "FUME"
    encoded_word_bytes = base64.b64encode(word_to_encode.encode('utf-8'))
    encoded_word = encoded_word_bytes.decode('utf-8')

    final_vector = np.array([0, 0, 0])
    vector_sum_parts = []
    
    for char in encoded_word:
        if 'A' <= char <= 'Z':
            vector_tuple = items[char]
            final_vector += np.array(vector_tuple)
            vector_sum_parts.append(str(vector_tuple))
            
    equation_str = " + ".join(vector_sum_parts)
    result_str = f"({final_vector[0]},{final_vector[1]},{final_vector[2]})"
    
    print(f"{equation_str} = {result_str}")

if __name__ == "__main__":
    solve_recommender_puzzle()
<<<(1,2,3)>>>