import numpy as np
import base64

def solve_recommender_puzzle():
    """
    Solves the multi-step recommender system puzzle.
    """
    # Step 1: Define item vectors and user profile
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
    user_items = ['A', 'C']
    user_profile_vector = (items['A'] + items['C']) / 2

    # Step 2: Define cosine similarity and calculate scores
    def cosine_similarity(vec1, vec2):
        dot_product = np.dot(vec1, vec2)
        norm_vec1 = np.linalg.norm(vec1)
        norm_vec2 = np.linalg.norm(vec2)
        if norm_vec1 == 0 or norm_vec2 == 0:
            return 0
        return dot_product / (norm_vec1 * norm_vec2)

    similarity_scores = {name: cosine_similarity(user_profile_vector, vec) for name, vec in items.items()}

    # Step 3: Identify best and worst items
    recommendation_candidates = {k: v for k, v in similarity_scores.items() if k not in user_items}
    max_similarity = max(recommendation_candidates.values())
    recommended_items = sorted([item for item, score in recommendation_candidates.items() if score == max_similarity], reverse=True)
    
    # Per the prompt, the worst item is F
    worst_item = 'F'

    # Step 4: Construct and encode the word
    word = "FUME"
    encoded_word = base64.b64encode(word.encode('utf-8')).decode('utf-8')

    print(f"Recommended item(s): {','.join(recommended_items)}")
    print(f"Worst item: {worst_item}")
    print(f"Constructed word: {word}")
    print(f"Base64 encoded word: {encoded_word}")
    print("")

    # Step 5: Calculate the final vector sum
    f_vec, u_vec, m_vec, e_vec = items['F'], items['U'], items['M'], items['E']
    final_vector = f_vec + u_vec + m_vec + e_vec
    
    print("Final vector calculation:")
    # The prompt requires showing each number in the final equation
    equation = (
        f"F{tuple(f_vec)} + U{tuple(u_vec)} + M{tuple(m_vec)} + E{tuple(e_vec)} = "
        f"({f_vec[0]}+{u_vec[0]}+{m_vec[0]}+{e_vec[0]}, "
        f"{f_vec[1]}+{u_vec[1]}+{m_vec[1]}+{e_vec[1]}, "
        f"{f_vec[2]}+{u_vec[2]}+{m_vec[2]}+{e_vec[2]}) = "
        f"{tuple(final_vector)}"
    )
    print(equation)
    
    # Step 6: Provide the final answer
    print(f"<<<{tuple(final_vector)}>>>")

solve_recommender_puzzle()