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

    # Step 1: Calculate the user's average feature vector
    user_interacted_items = ['A', 'C']
    user_vector_sum = np.sum([items[item] for item in user_interacted_items], axis=0)
    user_profile_vector = user_vector_sum / len(user_interacted_items)

    # Step 2: Calculate cosine similarity for all other items
    candidate_items = sorted([item for item in items if item not in user_interacted_items])
    similarities = {}

    norm_user = np.linalg.norm(user_profile_vector)

    for item_name in candidate_items:
        item_vector = items[item_name]
        norm_item = np.linalg.norm(item_vector)
        
        if norm_user == 0 or norm_item == 0:
            similarity = 0.0
        else:
            dot_product = np.dot(user_profile_vector, item_vector)
            similarity = dot_product / (norm_user * norm_item)
        
        similarities[item_name] = similarity

    # Step 3: Find the best and worst items
    max_similarity = -1
    if similarities:
        max_similarity = max(similarities.values())
    
    min_similarity = float('inf')
    if similarities:
        min_similarity = min(similarities.values())

    best_items = [name for name, sim in similarities.items() if sim == max_similarity]
    worst_items = [name for name, sim in similarities.items() if sim == min_similarity]

    # Apply tie-breaking rules
    best_items_sorted = sorted(best_items, reverse=True)
    
    worst_item_final = ''
    if 'F' in worst_items:
        worst_item_final = 'F'
    else:
        # Fallback in case 'F' isn't the worst, though it is in this problem
        worst_item_final = sorted(worst_items)[0] 

    # Step 4: Construct and encode the word
    # The letters are the best recommendations (U, M, E) and the worst (F).
    # A common English word from these letters is "FUME".
    word_letters = best_items_sorted + [worst_item_final]
    word = "FUME"
    encoded_word = base64.b64encode(word.encode('utf-8')).decode('utf-8')

    # Step 5: Sum the vectors of the letters in the word "FUME"
    final_vector_sum = np.sum([items[letter] for letter in word], axis=0)

    # Step 6: Print the results
    print(f"Recommended Items: {','.join(best_items_sorted)}")
    print(f"Worst Item: {worst_item_final}")
    print(f"Constructed Word: {word}")
    print(f"Base64 Encoded Word: {encoded_word}")

    sum_equation_parts = []
    for letter in word:
        vec_str = ','.join(map(str, items[letter]))
        sum_equation_parts.append(f"{letter}({vec_str})")
    
    final_vec_str = ','.join(map(str, final_vector_sum))
    print(f"Vector Summation: {' + '.join(sum_equation_parts)} = ({final_vec_str})")
    
    # Final Answer
    final_answer = f"({','.join(map(str, final_vector_sum))})"
    print(f"\n<<<{final_answer}>>>")


solve_recommender_problem()