import numpy as np
import base64
from collections import defaultdict

def solve_recommender_puzzle():
    """
    Solves the multi-step recommender system puzzle.
    """
    # Step 0: Define the item vectors
    items = {
        'A': (1, 0, 1), 'B': (0, 1, 1), 'C': (1, 1, 0), 'D': (0, 0, 1),
        'E': (1, 1, 1), 'F': (0, 0, 0), 'G': (1, 0, 0), 'H': (0, 1, 0),
        'I': (1, 0, 1), 'J': (0, 1, 1), 'K': (1, 1, 0), 'L': (0, 0, 1),
        'M': (1, 1, 1), 'N': (0, 0, 0), 'O': (1, 0, 0), 'P': (0, 1, 0),
        'Q': (1, 0, 1), 'R': (0, 1, 1), 'S': (1, 1, 0), 'T': (0, 0, 1),
        'U': (1, 1, 1), 'V': (0, 0, 0), 'W': (1, 0, 0), 'X': (0, 1, 0),
        'Y': (1, 0, 1), 'Z': (0, 1, 1)
    }

    # Convert to numpy arrays for easier calculation
    item_vectors = {name: np.array(vec) for name, vec in items.items()}

    # Step 1: Calculate the user's profile vector
    user_interacted_items = ['A', 'C']
    vec_a = item_vectors['A']
    vec_c = item_vectors['C']
    user_profile_vector = (vec_a + vec_c) / 2.0
    print(f"Step 1: The user's profile vector is the average of A={vec_a} and C={vec_c}, which is {user_profile_vector}\n")

    # Step 2: Calculate Cosine Similarities
    print("Step 2: Calculating cosine similarity with other items...")
    candidate_items = sorted([item for item in item_vectors if item not in user_interacted_items])
    
    similarities = {}
    
    # Helper for cosine similarity
    def cosine_similarity(v1, v2):
        norm_v1 = np.linalg.norm(v1)
        norm_v2 = np.linalg.norm(v2)
        if norm_v1 == 0 or norm_v2 == 0:
            return 0.0
        return np.dot(v1, v2) / (norm_v1 * norm_v2)

    for item in candidate_items:
        similarities[item] = cosine_similarity(user_profile_vector, item_vectors[item])
        # print(f"  - Similarity with {item} {item_vectors[item]}: {similarities[item]:.4f}")

    # Step 3: Find the best and worst items
    max_similarity = -1.0
    min_similarity = 2.0
    for sim in similarities.values():
        if sim > max_similarity:
            max_similarity = sim
        if sim < min_similarity:
            min_similarity = sim
            
    recommended_items = sorted([item for item, sim in similarities.items() if np.isclose(sim, max_similarity)], reverse=True)
    worst_items = [item for item, sim in similarities.items() if np.isclose(sim, min_similarity)]

    # Choose the worst item as per the prompt's cryptic instruction
    worst_item = 'F'
    
    print(f"\nStep 3: Identification of best and worst items.")
    print(f"  - Recommended item(s) (highest similarity): {','.join(recommended_items)}")
    print(f"  - Worst item(s) (lowest similarity): {', '.join(worst_items)}")
    print(f"  - Chosen worst item (least desired grade): {worst_item}")
    
    final_letters_string = f"{','.join(recommended_items)} {worst_item}"
    print(f"\nFinal letters for the next step: {final_letters_string}\n")
    
    # Step 4: Form word and encode in Base64
    word_to_encode = "FUME" # Anagram of E, F, M, U
    word_bytes = word_to_encode.encode('ascii')
    base64_encoded = base64.b64encode(word_bytes).decode('ascii')
    
    print(f"Step 4: Forming an English word from {{U,M,E,F}} and encoding it.")
    print(f"  - The word is '{word_to_encode}'.")
    print(f"  - Base64 encoding of '{word_to_encode}' is '{base64_encoded}'.\n")

    # Step 5: Sum the vectors from the Base64 result
    print("Step 5: Summing the vectors for letters from the Base64 string.")
    base64_letters = base64_encoded.replace('=', '')
    sum_vector = np.array([0, 0, 0])
    
    sum_x_str = []
    sum_y_str = []
    sum_z_str = []
    
    for char in base64_letters:
        # Treat lowercase as uppercase for vector lookup
        item_key = char.upper()
        vector = item_vectors[item_key]
        print(f"  - Letter '{char}' (maps to '{item_key}') has vector: {tuple(vector)}")
        sum_vector += vector
        sum_x_str.append(str(vector[0]))
        sum_y_str.append(str(vector[1]))
        sum_z_str.append(str(vector[2]))

    print("\nFinal equation:")
    print(f"X component = {' + '.join(sum_x_str)} = {sum_vector[0]}")
    print(f"Y component = {' + '.join(sum_y_str)} = {sum_vector[1]}")
    print(f"Z component = {' + '.join(sum_z_str)} = {sum_vector[2]}")

    final_vector = tuple(sum_vector)
    print(f"\nThe resulting final vector is {final_vector}.")
    
    # Final answer in the required format
    print(f"<<<{final_vector}>>>")

solve_recommender_puzzle()