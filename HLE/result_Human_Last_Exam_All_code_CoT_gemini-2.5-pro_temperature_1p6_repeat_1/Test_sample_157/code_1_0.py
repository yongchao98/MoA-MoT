import numpy as np
import base64
import re

def solve():
    """
    Solves the multi-step recommender system problem.
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

    # Step 1: Compute user's average feature vector
    user_interacted_items = ['A', 'C']
    user_vector = np.mean([items[item] for item in user_interacted_items], axis=0)
    print(f"User's average feature vector (profile) from items {user_interacted_items}: {tuple(np.round(user_vector, 2))}\n")


    # Step 2: Calculate cosine similarity for all other items
    def cosine_similarity(vec1, vec2):
        dot_product = np.dot(vec1, vec2)
        norm_vec1 = np.linalg.norm(vec1)
        norm_vec2 = np.linalg.norm(vec2)
        if norm_vec1 == 0 or norm_vec2 == 0:
            return 0.0
        return dot_product / (norm_vec1 * norm_vec2)

    candidate_items = sorted(list(set(items.keys()) - set(user_interacted_items)))
    
    similarities = {}
    for item in candidate_items:
        similarities[item] = cosine_similarity(user_vector, items[item])

    # Step 3: Find the item(s) with the highest and lowest cosine similarity
    max_similarity = -1.0
    for item in similarities:
        if similarities[item] > max_similarity:
            max_similarity = similarities[item]

    min_similarity = 2.0
    for item in similarities:
        if similarities[item] < min_similarity:
            min_similarity = similarities[item]
            
    recommended_items = sorted([item for item, sim in similarities.items() if sim == max_similarity], reverse=True)
    worst_items_candidates = [item for item, sim in similarities.items() if sim == min_similarity]
    
    # Per prompt, choose 'F' if available among the worst options
    worst_item = 'F' if 'F' in worst_items_candidates else worst_items_candidates[0]
    
    print(f"Recommended item(s) (highest similarity): {','.join(recommended_items)}")
    print(f"Worst item (based on prompt rule): {worst_item}\n")

    # Step 4: Construct and encode the word
    best_letters_str = "".join(recommended_items)
    word_str = worst_item + best_letters_str[2] + best_letters_str[1] + best_letters_str[0] # FUME
    word_bytes = word_str.encode('ascii')
    encoded_word = base64.b64encode(word_bytes).decode('ascii')
    
    print(f"Constructed word: {word_str}")
    print(f"Base64 encoded word: {encoded_word}\n")

    # Step 5: Extract letters and sum their vectors
    letters_from_base64 = re.findall('[a-zA-Z]', encoded_word)
    # Treat lowercase as uppercase for vector lookup
    letters_to_sum = [char.upper() for char in letters_from_base64]
    
    print(f"Letters from Base64 string for vector sum: {', '.join(letters_to_sum)}")
    
    sum_vector = np.array([0, 0, 0])
    equation_parts = []
    for letter in letters_to_sum:
        vector = items[letter]
        sum_vector = sum_vector + vector
        equation_parts.append(f"{letter}{tuple(vector)}")
        
    final_equation = " + ".join(equation_parts)
    final_vector_tuple = tuple(sum_vector)

    print("Sum of vectors:")
    print(f"{final_equation} = {final_vector_tuple}")
    
    # Final answer
    print(f"\nFinal resulting vector: {final_vector_tuple}")
    return final_vector_tuple

if __name__ == "__main__":
    final_answer = solve()
    print(f"<<<{final_answer}>>>")