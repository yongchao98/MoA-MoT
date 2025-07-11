import numpy as np

def solve_recommender_problem():
    """
    Solves the recommender system problem by finding the best/worst items,
    forming a word, and summing their vectors.
    """
    # Define the item feature vectors
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
    
    # Items the user has interacted with
    user_items = ['A', 'C']

    # Step 1: Compute the average feature vector for the user
    user_vectors = [items[item] for item in user_items]
    user_profile_vector = np.mean(user_vectors, axis=0)

    def cosine_similarity(v1, v2):
        norm_v1 = np.linalg.norm(v1)
        norm_v2 = np.linalg.norm(v2)
        if norm_v1 == 0 or norm_v2 == 0:
            return 0.0
        return np.dot(v1, v2) / (norm_v1 * norm_v2)

    # Step 2: Calculate similarities for candidate items
    similarities = {}
    recommendation_candidates = {k: v for k, v in items.items() if k not in user_items}

    for item, vector in recommendation_candidates.items():
        sim = cosine_similarity(user_profile_vector, vector)
        similarities[item] = sim

    # Step 3: Find the best and worst items
    max_similarity = max(similarities.values())
    min_similarity = min(similarities.values())

    best_items = sorted([item for item, sim in similarities.items() if sim == max_similarity], reverse=True)
    worst_items = [item for item, sim in similarities.items() if sim == min_similarity]
    
    # Per instructions, the worst letter is 'F'
    worst_letter = 'F'
    
    # The letters for the word are the best items (U, M, E) and the worst item (F)
    # Constructing the word "FUME"
    word_letters = [worst_letter] + best_items

    # Step 4: Add the vectors of the letters in the word "FUME"
    vectors_to_sum = [items[letter] for letter in word_letters]
    result_vector = np.sum(vectors_to_sum, axis=0)
    
    # Step 5: Output the final equation
    equation_parts = []
    for vec in vectors_to_sum:
        equation_parts.append(f"({vec[0]}, {vec[1]}, {vec[2]})")
    
    equation_str = " + ".join(equation_parts)
    result_str = f"({result_vector[0]}, {result_vector[1]}, {result_vector[2]})"
    
    print(f"{equation_str} = {result_str}")

solve_recommender_problem()