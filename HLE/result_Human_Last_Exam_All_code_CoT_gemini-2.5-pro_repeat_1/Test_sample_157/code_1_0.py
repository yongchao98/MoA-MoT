import numpy as np
import base64

def solve_recommender_task():
    """
    Solves the multi-step recommender system problem.
    """
    # 1. Define item vectors
    items = {
        'A': (1, 0, 1), 'B': (0, 1, 1), 'C': (1, 1, 0), 'D': (0, 0, 1),
        'E': (1, 1, 1), 'F': (0, 0, 0), 'G': (1, 0, 0), 'H': (0, 1, 0),
        'I': (1, 0, 1), 'J': (0, 1, 1), 'K': (1, 1, 0), 'L': (0, 0, 1),
        'M': (1, 1, 1), 'N': (0, 0, 0), 'O': (1, 0, 0), 'P': (0, 1, 0),
        'Q': (1, 0, 1), 'R': (0, 1, 1), 'S': (1, 1, 0), 'T': (0, 0, 1),
        'U': (1, 1, 1), 'V': (0, 0, 0), 'W': (1, 0, 0), 'X': (0, 1, 0),
        'Y': (1, 0, 1), 'Z': (0, 1, 1)
    }
    item_vectors = {key: np.array(val) for key, val in items.items()}

    # 2. Calculate the user's profile vector
    user_items = ['A', 'C']
    user_vectors = [item_vectors[item] for item in user_items]
    profile_vector = np.mean(user_vectors, axis=0)
    print(f"Step 1: The user's profile vector (average of A and C) is {tuple(np.round(profile_vector, 2))}")

    # 3. Compute cosine similarities for candidate items
    def cosine_similarity(v1, v2):
        dot_product = np.dot(v1, v2)
        norm_v1 = np.linalg.norm(v1)
        norm_v2 = np.linalg.norm(v2)
        if norm_v1 == 0 or norm_v2 == 0:
            return 0.0
        return dot_product / (norm_v1 * norm_v2)

    candidate_items = sorted([item for item in item_vectors.keys() if item not in user_items])
    similarities = {item: cosine_similarity(profile_vector, item_vectors[item]) for item in candidate_items}
    print("\nStep 2: Cosine similarities with the user profile have been calculated for all other items.")

    # 4. Identify best and worst items
    max_sim = max(similarities.values())
    min_sim = min(similarities.values())
    
    best_items = sorted([item for item, sim in similarities.items() if sim == max_sim], reverse=True)
    # Per instructions, choose 'F' from the items with the lowest similarity
    worst_item = 'F'
    
    print(f"\nStep 3: Identification of key items.")
    print(f"  - Recommended items (highest similarity, reverse lexicographical): {', '.join(best_items)}")
    print(f"  - Worst item (lowest similarity, specific choice): {worst_item}")
    
    # 5. Construct and encode the word
    # The letters are from the best items (U, M, E) and the worst item (F).
    # These letters can form the English word "FUME".
    word = "FUME"
    encoded_word = base64.b64encode(word.encode('ascii')).decode('ascii')
    print(f"\nStep 4: The letters U, M, E, F form the word '{word}'.")
    print(f"  - The Base64 encoding of '{word}' is: {encoded_word}")

    # 6. Calculate the final vector sum
    sum_vector = np.array([0, 0, 0])
    equation_parts = []
    for letter in word:
        vec = item_vectors[letter]
        sum_vector += vec
        
    print(f"\nStep 5: The vectors for the letters in '{word}' are summed up.")
    f_vec = item_vectors['F']
    u_vec = item_vectors['U']
    m_vec = item_vectors['M']
    e_vec = item_vectors['E']
    s_vec = sum_vector
    
    # Final equation printing each number
    print(f"The final equation is: Vector F({f_vec[0]}, {f_vec[1]}, {f_vec[2]}) + Vector U({u_vec[0]}, {u_vec[1]}, {u_vec[2]}) + Vector M({m_vec[0]}, {m_vec[1]}, {m_vec[2]}) + Vector E({e_vec[0]}, {e_vec[1]}, {e_vec[2]}) = ({s_vec[0]}, {s_vec[1]}, {s_vec[2]})")

    # 7. Final Answer
    print(f"\n<<<{encoded_word}>>>")

if __name__ == '__main__':
    solve_recommender_task()