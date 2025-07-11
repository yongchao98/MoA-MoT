import numpy as np
import base64

def solve_recommender_problem():
    """
    This function solves the multi-step recommender system problem.
    """
    # Step 1: Define the item vectors
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

    # Step 2: Calculate the user's average profile vector
    user_interacted_items = ['A', 'C']
    user_profile_vector = np.mean([items[item] for item in user_interacted_items], axis=0)

    # Define a helper function for cosine similarity
    def cosine_similarity(v1, v2):
        dot_product = np.dot(v1, v2)
        norm_v1 = np.linalg.norm(v1)
        norm_v2 = np.linalg.norm(v2)
        # Avoid division by zero for zero vectors
        if norm_v1 == 0 or norm_v2 == 0:
            return 0.0
        return dot_product / (norm_v1 * norm_v2)

    # Step 3: Calculate similarities for all candidate items (not previously interacted with)
    candidate_items = {key: val for key, val in items.items() if key not in user_interacted_items}
    similarities = {item: cosine_similarity(user_profile_vector, vec) for item, vec in candidate_items.items()}

    # Find the maximum and minimum similarity scores
    max_similarity = max(similarities.values())
    min_similarity = min(similarities.values())

    # Identify all items with the highest similarity
    recommended_items = sorted([item for item, sim in similarities.items() if np.isclose(sim, max_similarity)], reverse=True)
    
    # Identify all items with the lowest similarity
    worst_items = [item for item, sim in similarities.items() if np.isclose(sim, min_similarity)]
    # Select 'F' as the worst item per the instructions
    worst_item_letter = 'F'

    # Step 4: Construct the word from the identified letters
    # The best items are E, M, U. The worst is F. A plausible word is "FUME".
    word_letters = ['F', 'U', 'M', 'E']
    word = "".join(word_letters)
    
    # Step 5: Sum the vectors for the letters in the word "FUME"
    sum_vector = np.zeros(3, dtype=int)
    equation_parts = []
    for letter in word_letters:
        vec = items[letter]
        sum_vector += vec
        equation_parts.append(f"{letter}{tuple(vec)}")

    # Step 6: Print the final equation showing each number
    print("The final calculation is the sum of the vectors for the letters in the word 'FUME':")
    
    x_nums = "+".join(str(items[l][0]) for l in word_letters)
    y_nums = "+".join(str(items[l][1]) for l in word_letters)
    z_nums = "+".join(str(items[l][2]) for l in word_letters)

    equation_str = f"{' + '.join(equation_parts)} = ({x_nums}, {y_nums}, {z_nums}) = {tuple(sum_vector)}"
    print(equation_str)

    # Output the final answer in the specified format
    final_answer = tuple(sum_vector)
    print(f"<<<{final_answer}>>>")

# Execute the function
solve_recommender_problem()