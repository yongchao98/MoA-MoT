import numpy as np
import base64

def solve_recommender_problem():
    """
    Solves the multi-step recommender system problem.
    """
    items_data = {
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

    # 1. Calculate the user's average profile vector
    user_items = ['A', 'C']
    user_vectors = [items_data[item] for item in user_items]
    user_profile_vector = np.mean(user_vectors, axis=0)
    
    # 2. Define the cosine similarity function
    def cosine_similarity(vec1, vec2):
        dot_product = np.dot(vec1, vec2)
        norm_vec1 = np.linalg.norm(vec1)
        norm_vec2 = np.linalg.norm(vec2)
        if norm_vec1 == 0 or norm_vec2 == 0:
            return 0.0
        return dot_product / (norm_vec1 * norm_vec2)

    # 3. Calculate similarities for all candidate items
    candidate_items = [item for item in items_data if item not in user_items]
    similarities = {
        item: cosine_similarity(user_profile_vector, items_data[item])
        for item in candidate_items
    }

    # 4. Find the best and worst items
    if not similarities:
        print("No candidate items to recommend.")
        return

    max_similarity = max(similarities.values())
    min_similarity = min(similarities.values())

    best_items = sorted(
        [item for item, sim in similarities.items() if np.isclose(sim, max_similarity)],
        reverse=True
    )
    
    # Per instructions, choose 'F' from the worst items.
    worst_item = 'F'

    print(f"Recommended item(s): {', '.join(best_items)}")
    print(f"Worst item: {worst_item}")
    
    # 5. Form a word, encode it, and sum vectors
    word_letters = best_items + [worst_item]
    # Rearrange to form a recognizable English word
    word = "FUME"
    encoded_word = base64.b64encode(word.encode('utf-8')).decode('utf-8')
    
    print(f"Formed word: {word}")
    print(f"Base64 encoded word: {encoded_word}")

    sum_vector = np.sum([items_data[letter] for letter in word], axis=0)

    f_vec = tuple(items_data['F'])
    u_vec = tuple(items_data['U'])
    m_vec = tuple(items_data['M'])
    e_vec = tuple(items_data['E'])
    
    print(f"Vector sum equation: F{f_vec} + U{u_vec} + M{m_vec} + E{e_vec} = {tuple(sum_vector)}")
    
    # Final answer
    print(f"<<<{encoded_word}>>>")


solve_recommender_problem()