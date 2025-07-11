import numpy as np
import base64

def solve_recommender_task():
    """
    Solves the recommender system problem by calculating the best and worst recommendations,
    forming a word, and summing the corresponding vectors.
    """
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

    # Step 1: Calculate the average user profile vector
    user_items = ['A', 'C']
    user_vectors = [items[item] for item in user_items]
    user_profile_sum = np.sum(user_vectors, axis=0)
    user_profile_avg = user_profile_sum / len(user_items)
    
    print("--- User Profile Calculation ---")
    print(f"User interacted with: {user_items}")
    print(f"Vector for A: {items['A']}")
    print(f"Vector for C: {items['C']}")
    print(f"Sum of vectors: {items['A']} + {items['C']} = {user_profile_sum}")
    print(f"Average user profile vector: {user_profile_sum} / {len(user_items)} = {user_profile_avg}\n")

    # Step 2: Define candidate items and calculate cosine similarity
    candidate_items = sorted([item for item in items if item not in user_items])
    
    def cosine_similarity(vec1, vec2):
        dot_product = np.dot(vec1, vec2)
        norm_vec1 = np.linalg.norm(vec1)
        norm_vec2 = np.linalg.norm(vec2)
        if norm_vec1 == 0 or norm_vec2 == 0:
            return 0.0 # Similarity with a zero vector is 0
        return dot_product / (norm_vec1 * norm_vec2)

    similarities = {}
    print("--- Cosine Similarity Scores with User Profile ---")
    for item in candidate_items:
        sim = cosine_similarity(user_profile_avg, items[item])
        similarities[item] = sim
        print(f"Similarity with '{item}' {items[item]}: {sim:.4f}")
    
    # Step 3: Find the best and worst items
    max_sim = max(similarities.values())
    min_sim = min(similarities.values())
    
    best_items = [item for item, sim in similarities.items() if sim == max_sim]
    best_items.sort(reverse=True) # Sort in reverse lexicographical order

    worst_items = [item for item, sim in similarities.items() if sim == min_sim]
    # As per the prompt, choose 'F' for the worst item.
    worst_item = 'F'
    
    print("\n--- Recommendation Results ---")
    print(f"Highest similarity score: {max_sim:.4f}")
    print(f"Recommended items (in reverse lexicographical order): {', '.join(best_items)}")
    print(f"Lowest similarity score: {min_sim:.4f}")
    print(f"Worst item (based on rule): {worst_item}\n")

    # Step 4: Form a word and encode it
    word_letters = best_items + [worst_item]
    # Anagram of {'U', 'M', 'E', 'F'} is 'FUME'
    word = "FUME"
    word_uppercase = word.upper()
    encoded_word = base64.b64encode(word_uppercase.encode('utf-8'))
    
    print("--- Word Construction and Encoding ---")
    print(f"Letters for word construction: {word_letters}")
    print(f"Constructed English word: {word_uppercase}")
    print(f"Base64 encoded word: {encoded_word.decode('utf-8')}\n")

    # Step 5: Sum the vectors of the letters in the constructed word
    final_vector_sum = np.array([0, 0, 0])
    vector_components = []
    for letter in word_uppercase:
        final_vector_sum += items[letter]
        vector_components.append(str(items[letter]))

    print("--- Final Vector Calculation ---")
    equation = " + ".join(vector_components)
    print(f"Final Equation: {equation} = {final_vector_sum}")

    # Final answer format
    final_answer = f"({final_vector_sum[0]}, {final_vector_sum[1]}, {final_vector_sum[2]})"
    print(f"\n<<<{final_answer}>>>")

solve_recommender_task()