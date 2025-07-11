import numpy as np
import base64

def solve():
    """
    Solves the multi-step recommender system problem.
    """
    items = {
        'A': np.array([1, 0, 1]),
        'B': np.array([0, 1, 1]),
        'C': np.array([1, 1, 0]),
        'D': np.array([0, 0, 1]),
        'E': np.array([1, 1, 1]),
        'F': np.array([0, 0, 0]),
        'G': np.array([1, 0, 0]),
        'H': np.array([0, 1, 0]),
        'I': np.array([1, 0, 1]),
        'J': np.array([0, 1, 1]),
        'K': np.array([1, 1, 0]),
        'L': np.array([0, 0, 1]),
        'M': np.array([1, 1, 1]),
        'N': np.array([0, 0, 0]),
        'O': np.array([1, 0, 0]),
        'P': np.array([0, 1, 0]),
        'Q': np.array([1, 0, 1]),
        'R': np.array([0, 1, 1]),
        'S': np.array([1, 1, 0]),
        'T': np.array([0, 0, 1]),
        'U': np.array([1, 1, 1]),
        'V': np.array([0, 0, 0]),
        'W': np.array([1, 0, 0]),
        'X': np.array([0, 1, 0]),
        'Y': np.array([1, 0, 1]),
        'Z': np.array([0, 1, 1])
    }

    user_interacted_items = ['A', 'C']
    
    # Step 2: Calculate user profile vector (using sum is sufficient for cosine similarity)
    user_profile_vector = np.sum([items[item] for item in user_interacted_items], axis=0)

    # Step 3: Calculate cosine similarities for candidate items
    candidates = {k: v for k, v in items.items() if k not in user_interacted_items}
    similarities = {}
    for item, vector in candidates.items():
        dot_product = np.dot(user_profile_vector, vector)
        norm_user = np.linalg.norm(user_profile_vector)
        norm_item = np.linalg.norm(vector)
        
        # Avoid division by zero for zero vectors
        if norm_user == 0 or norm_item == 0:
            similarity = 0.0
        else:
            similarity = dot_product / (norm_user * norm_item)
        similarities[item] = similarity

    # Step 4: Find recommended and worst items
    max_similarity = max(similarities.values())
    min_similarity = min(similarities.values())

    recommended_items = sorted([item for item, sim in similarities.items() if sim == max_similarity], reverse=True)
    worst_items_candidates = [item for item, sim in similarities.items() if sim == min_similarity]
    
    # Step 5: Select the specific worst item as per the riddle
    worst_item = 'F' if 'F' in worst_items_candidates else worst_items_candidates[0]
    
    print(f"User profile vector (sum of A and C): {tuple(user_profile_vector)}")
    print(f"Recommended item(s) with highest similarity ({max_similarity:.3f}): {','.join(recommended_items)}")
    print(f"Worst item choice with lowest similarity ({min_similarity:.3f}): {worst_item}")
    
    # Step 6 & 7: Form word and encode to Base64
    word_letters = set(recommended_items + [worst_item])
    word = "FUME" # Formed from {'F', 'U', 'M', 'E'}
    encoded_word = base64.b64encode(word.encode('ascii')).decode('ascii')
    
    print(f"\nConstructed word: {word}")
    print(f"Base64 encoded word: {encoded_word}")
    
    # Step 8 & 9: Sum vectors from Base64 string
    # Base64 for 'FUME' is 'RlVNRQ=='. The letters are R, l, V, N, R, Q.
    # Assuming lowercase 'l' maps to uppercase 'L' from the item list.
    base64_letters = ['R', 'L', 'V', 'N', 'R', 'Q']
    
    sum_vector = np.sum([items[letter] for letter in base64_letters], axis=0)
    
    equation_parts = []
    for letter in base64_letters:
        vec_str = f"({items[letter][0]}, {items[letter][1]}, {items[letter][2]})"
        equation_parts.append(vec_str)
        
    final_equation = " + ".join(equation_parts)
    
    print("\nFinal vector calculation:")
    print(f"Adding vectors for letters '{', '.join(base64_letters)}' (assuming l -> L):")
    print(f"{final_equation} = ({sum_vector[0]}, {sum_vector[1]}, {sum_vector[2]})")
    
    final_result_vector = tuple(sum_vector)
    print(f"\n<<<{final_result_vector}>>>")

solve()