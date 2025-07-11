import numpy as np
import base64

def solve_recommender_problem():
    """
    Solves the multi-step recommender system problem.
    """
    # Step 1: Define item vectors
    items = {
        'A': (1, 0, 1), 'B': (0, 1, 1), 'C': (1, 1, 0), 'D': (0, 0, 1),
        'E': (1, 1, 1), 'F': (0, 0, 0), 'G': (1, 0, 0), 'H': (0, 1, 0),
        'I': (1, 0, 1), 'J': (0, 1, 1), 'K': (1, 1, 0), 'L': (0, 0, 1),
        'M': (1, 1, 1), 'N': (0, 0, 0), 'O': (1, 0, 0), 'P': (0, 1, 0),
        'Q': (1, 0, 1), 'R': (0, 1, 1), 'S': (1, 1, 0), 'T': (0, 0, 1),
        'U': (1, 1, 1), 'V': (0, 0, 0), 'W': (1, 0, 0), 'X': (0, 1, 0),
        'Y': (1, 0, 1), 'Z': (0, 1, 1)
    }

    # Step 2: Define user interactions and compute user profile vector
    interacted_items = ['A', 'C']
    interacted_vectors = [np.array(items[item]) for item in interacted_items]
    user_profile_vector = np.mean(interacted_vectors, axis=0)

    # Cosine similarity function
    def cosine_similarity(v1, v2):
        dot_product = np.dot(v1, v2)
        norm_v1 = np.linalg.norm(v1)
        norm_v2 = np.linalg.norm(v2)
        if norm_v1 == 0 or norm_v2 == 0:
            return 0
        return dot_product / (norm_v1 * norm_v2)

    # Step 3: Calculate similarities for candidate items
    candidate_items = {k: v for k, v in items.items() if k not in interacted_items}
    similarities = {
        item: cosine_similarity(user_profile_vector, np.array(vector))
        for item, vector in candidate_items.items()
    }

    # Step 4: Find the best and worst items
    if not similarities:
        print("No candidate items to recommend.")
        return

    max_similarity = max(similarities.values())
    min_similarity = min(similarities.values())

    best_items = sorted([item for item, sim in similarities.items() if sim == max_similarity], reverse=True)
    worst_items_options = {item for item, sim in similarities.items() if sim == min_similarity}
    
    # Choose the worst item based on the prompt's clue
    worst_item = 'F' if 'F' in worst_items_options else sorted(list(worst_items_options))[0]

    # Step 5: Form a word and encode it
    # The letters are E, M, U from best_items and F from worst_item.
    # A common English word from {E, F, M, U} is FUME.
    word = "FUME"
    encoded_word = base64.b64encode(word.encode('ascii')).decode('ascii')
    
    # Step 6: Sum vectors from Base64 string letters
    base64_letters = [char for char in encoded_word if char.isalpha()]
    
    sum_vector = np.array([0, 0, 0])
    vector_list_str = []
    
    for letter in base64_letters:
        # Treat lowercase as uppercase for vector lookup
        vector = np.array(items[letter.upper()])
        sum_vector += vector
        vector_list_str.append(str(tuple(vector)))

    # Step 7: Print the final equation and the answer
    equation = " + ".join(vector_list_str)
    final_vector_tuple = tuple(sum_vector)

    print(f"The best recommended items are: {', '.join(best_items)}")
    print(f"The worst item is: {worst_item}")
    print(f"The constructed word is: {word}")
    print(f"The Base64 encoding is: {encoded_word}")
    print(f"The letters from the Base64 string are: {', '.join(base64_letters)}")
    print("\nCalculating the final vector sum:")
    print(f"{equation} = {final_vector_tuple}")
    
    # Final answer in the required format
    print(f"\n<<<{final_vector_tuple}>>>")

solve_recommender_problem()