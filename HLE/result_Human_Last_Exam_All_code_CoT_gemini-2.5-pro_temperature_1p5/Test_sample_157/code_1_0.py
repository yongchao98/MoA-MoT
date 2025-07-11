import numpy as np
import base64

def solve_recommender_puzzle():
    """
    Solves the entire multi-step recommender system problem.
    """
    # Step 1: Define data and calculate user profile vector
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

    user_interacted_items = ['A', 'C']
    vec_A = items['A']
    vec_C = items['C']

    user_profile_vector = (vec_A + vec_C) / 2.0
    print(f"User interacted with A={tuple(vec_A)} and C={tuple(vec_C)}")
    print(f"User profile vector (average) = {tuple(user_profile_vector)}\n")

    # Step 2: Define and compute Cosine Similarities
    def cosine_similarity(v1, v2):
        dot_product = np.dot(v1, v2)
        norm_v1 = np.linalg.norm(v1)
        norm_v2 = np.linalg.norm(v2)
        if norm_v1 == 0 or norm_v2 == 0:
            return 0.0
        return dot_product / (norm_v1 * norm_v2)

    candidate_items = {name: vector for name, vector in items.items() if name not in user_interacted_items}

    similarities = {name: cosine_similarity(user_profile_vector, vector) for name, vector in candidate_items.items()}

    # Step 3: Identify Best and Worst Items
    max_sim = max(similarities.values())
    min_sim = min(similarities.values())

    # Get all items with the max similarity and sort them in reverse lexicographical order
    best_items = sorted([name for name, sim in similarities.items() if np.isclose(sim, max_sim)], reverse=True)
    
    # Get all items with the min similarity and choose 'F' per the instructions
    worst_item_choice = 'F'
    
    print(f"Highest similarity score found: {max_sim:.4f}")
    print(f"Recommended items (in reverse lexicographical order): {', '.join(best_items)}")
    print(f"Chosen worst item based on rule: {worst_item_choice}\n")

    # Step 4: Form Word and Encode
    # The letters for the word are from the best (U, M, E) and worst (F) items.
    # The English word formed is FUME.
    word = "FUME"
    print(f"Constructed word: {word}")
    
    encoded_word = base64.b64encode(word.encode('ascii')).decode('ascii')
    print(f"Base64 encoded word: {encoded_word}\n")

    # Step 5: Calculate Final Vector by summing vectors of letters in the word "FUME"
    final_vector = np.array([0, 0, 0])
    equation_parts = []
    
    word_components = list(word) # ['F', 'U', 'M', 'E']
    
    print("Final vector calculation (sum of vectors for letters in FUME):")
    for letter in word_components:
        vector = items[letter]
        final_vector = final_vector + vector
        equation_parts.append(f"{letter}{tuple(vector)}")
        
    equation_str = " + ".join(equation_parts)
    print(f"{equation_str} = {tuple(final_vector)}")
    
    # Final Output
    print(f"\n<<<{tuple(final_vector)}>>>")

# Run the solution
solve_recommender_puzzle()