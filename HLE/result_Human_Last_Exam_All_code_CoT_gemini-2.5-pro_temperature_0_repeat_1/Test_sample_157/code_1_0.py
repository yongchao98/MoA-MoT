import numpy as np
import base64
from collections import defaultdict

def solve_recommender_puzzle():
    """
    Solves the multi-step recommender system puzzle.
    """
    # Step 1: Define the item vectors
    items = {
        'A': (1, 0, 1), 'B': (0, 1, 1), 'C': (1, 1, 0), 'D': (0, 0, 1),
        'E': (1, 1, 1), 'F': (0, 0, 0), 'G': (1, 0, 0), 'H': (0, 1, 0),
        'I': (1, 0, 1), 'J': (0, 1, 1), 'K': (1, 1, 0), 'L': (0, 0, 1),
        'M': (1, 1, 1), 'N': (0, 0, 0), 'O': (1, 0, 0), 'P': (0, 1, 0),
        'Q': (1, 0, 1), 'R': (0, 1, 1), 'S': (1, 1, 0), 'T': (0, 0, 1),
        'U': (1, 1, 1), 'V': (0, 0, 0), 'W': (1, 0, 0), 'X': (0, 1, 0),
        'Y': (1, 0, 1), 'Z': (0, 1, 1)
    }

    # Step 2: Calculate the user's profile vector
    user_items = ['A', 'C']
    user_vectors = [np.array(items[item]) for item in user_items]
    profile_vector = np.mean(user_vectors, axis=0)

    # Step 3: Calculate cosine similarity for all other items
    def cosine_similarity(v1, v2):
        dot_product = np.dot(v1, v2)
        norm_v1 = np.linalg.norm(v1)
        norm_v2 = np.linalg.norm(v2)
        # Handle zero vectors to avoid division by zero
        if norm_v1 == 0 or norm_v2 == 0:
            return 0.0
        return dot_product / (norm_v1 * norm_v2)

    similarities = {}
    candidate_items = sorted(list(items.keys()))
    for item_name in candidate_items:
        if item_name not in user_items:
            item_vector = np.array(items[item_name])
            sim = cosine_similarity(profile_vector, item_vector)
            similarities[item_name] = sim

    # Step 4: Find the best and worst recommendations
    # Group items by their similarity score to handle ties
    sim_groups = defaultdict(list)
    for item, sim in similarities.items():
        # Use a rounded similarity score as the key to group floats accurately
        sim_groups[round(sim, 8)].append(item)

    max_sim = max(sim_groups.keys())
    min_sim = min(sim_groups.keys())

    best_items = sim_groups[max_sim]
    best_items.sort(reverse=True)  # Sort in reverse lexicographical order

    worst_items_candidates = sim_groups[min_sim]
    # From {F, N, V}, 'F' is the grade students least want
    worst_item = 'F'

    # Step 5: Construct the word from the identified letters
    # The letters are from the best items {U, M, E} and the worst item {F}
    # Anagram to form an English word: FUME
    word = "FUME"

    # Step 6: Base64 encode the word
    encoded_word = base64.b64encode(word.encode('ascii')).decode('ascii')

    # Step 7: Sum the vectors for the letters in the word "FUME"
    word_vectors = [items[char] for char in word]
    sum_vector = np.sum(np.array(word_vectors), axis=0)

    # Step 8: Print all results and the final equation
    print(f"1. User profile vector (average of A and C): {tuple(profile_vector)}")
    print(f"2. Recommended item(s) (highest similarity): {', '.join(best_items)}")
    print(f"3. Worst item (lowest similarity, 'F' grade): {worst_item}")
    print(f"4. Constructed word: {word}")
    print(f"5. Base64 encoded word: {encoded_word}")
    
    print("\n6. Final vector sum calculation:")
    v_f, v_u, v_m, v_e = items['F'], items['U'], items['M'], items['E']
    
    # Print the equation with each number explicitly shown
    print(f"   Vector('F') + Vector('U') + Vector('M') + Vector('E')")
    print(f"=  {v_f} + {v_u} + {v_m} + {v_e}")
    print(f"=  ({v_f[0]}+{v_u[0]}+{v_m[0]}+{v_e[0]}, {v_f[1]}+{v_u[1]}+{v_m[1]}+{v_e[1]}, {v_f[2]}+{v_u[2]}+{v_m[2]}+{v_e[2]})")
    print(f"=  ({sum_vector[0]}, {sum_vector[1]}, {sum_vector[2]})")

# Execute the main function
solve_recommender_puzzle()