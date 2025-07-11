import numpy as np
import base64

def solve_recommender_problem():
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

    user_items = ['A', 'C']
    
    # Step 1: Calculate user's average feature vector
    user_vector_sum = np.sum([items[item] for item in user_items], axis=0)
    user_profile_vector = user_vector_sum / len(user_items)
    print(f"User Profile Vector (Average of A and C): {user_profile_vector}\n")

    # Step 2: Calculate cosine similarity with all other items
    items_to_recommend = sorted([item for item in items if item not in user_items])
    similarities = {}
    print("Cosine Similarities:")
    
    norm_user_profile = np.linalg.norm(user_profile_vector)

    for item in items_to_recommend:
        item_vector = items[item]
        norm_item = np.linalg.norm(item_vector)
        
        # Handle zero vector case
        if norm_item == 0 or norm_user_profile == 0:
            similarity = 0.0
        else:
            dot_product = np.dot(user_profile_vector, item_vector)
            similarity = dot_product / (norm_user_profile * norm_item)
        
        similarities[item] = similarity
        print(f"  - Sim(User, {item}) = {similarity:.4f}")

    # Step 3: Find the best and worst recommendations
    max_similarity = max(similarities.values())
    min_similarity = min(similarities.values())

    best_items = sorted([item for item, sim in similarities.items() if sim == max_similarity], reverse=True)
    worst_items_candidates = [item for item, sim in similarities.items() if sim == min_similarity]
    
    # Per instructions, choose 'F' if available among worst candidates
    worst_item = 'F' if 'F' in worst_items_candidates else worst_items_candidates[0]
    
    print("\n--- Results ---")
    print(f"Recommended Item(s) (Highest Similarity): {', '.join(best_items)}")
    print(f"Worst Item (Lowest Similarity): {worst_item}")
    
    # Step 4: Construct and encode the word
    word_letters = best_items + [worst_item]
    word = "".join(word_letters)
    word_uppercase = word.upper()
    
    encoded_word = base64.b64encode(word_uppercase.encode('utf-8')).decode('utf-8')
    
    print(f"Constructed Word: {word}")
    print(f"Base64 Encoded Word: {encoded_word}\n")
    
    # Step 5: Sum vectors of the letters in the word
    print("--- Final Vector Calculation ---")
    final_vector = np.array([0, 0, 0])
    
    equation_parts = []
    for letter in word_letters:
        vec = items[letter]
        final_vector += vec
        equation_parts.append(f"Vector({letter})")
        
    print(f"Sum Vector = {' + '.join(equation_parts)}")

    vector_strings = [str(tuple(items[letter])) for letter in word_letters]
    print(f"Sum Vector = {' + '.join(vector_strings)}")

    sum_x = "+".join([str(items[letter][0]) for letter in word_letters])
    sum_y = "+".join([str(items[letter][1]) for letter in word_letters])
    sum_z = "+".join([str(items[letter][2]) for letter in word_letters])

    print(f"Sum Vector = ({sum_x}, {sum_y}, {sum_z})")
    
    final_vector_tuple = tuple(final_vector)
    print(f"Final Resulting Vector = {final_vector_tuple}")

# Execute the function
solve_recommender_problem()